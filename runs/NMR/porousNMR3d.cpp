/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Fabian Klemens
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#define porousMode

#include "olb3D.h"
typedef double S;
typedef double T;
typedef double BaseType;

#include "olb3D.hh"
#include "controlled.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::opti;
using namespace olb::util;

#define Lattice DualPorousD3Q19Descriptor

// global analytical solution to optimize against
olb::AnalyticalF3D<T,T>* solution;

static unsigned int counterOptiStep = 0;
bool optimizationStarted = false;
static int __overlap = 2;
static std::vector<T> errorNorm;
static std::vector<T> errorNormInside;
static std::vector<int> numberOfSimulations;
static int errorNormCounter = 0;
std::string errorPlotFile = "errorPlot.dat";

SuperIndicatorF3D<T>* objectiveDomain;

string output_dir;
string getOptiDirName()
{
  return output_dir + "opti_primal_step_" + to_string(counterOptiStep) + "/";
}


template <typename T, template<typename U> class Lattice>
void error(SuperGeometry3D<T>& superGeometry,
           SuperLattice3D<T, Lattice>& lattice,
           UnitConverter<T,Lattice> const& converter)
{
  OstreamManager clout(std::cout,"error");

  T error_abs, u_abs;

  SuperLatticePhysVelocity3D<T,Lattice> u(lattice, converter);
  SuperLatticeFfromAnalyticalF3D<T,Lattice> uSolLattice(*solution, lattice);

  SuperL2Norm3D<T> diffL2Norm(uSolLattice-u, superGeometry, {1,5,6});
  SuperL2Norm3D<T> absL2Norm(uSolLattice, superGeometry, {1,5,6});
  diffL2Norm(&error_abs, nullptr);
  absL2Norm(&u_abs, nullptr);
  clout << "Mat. 1,5,6: ";
  clout << "velocity-L2-error(abs)=" << error_abs << "; \t";
  clout << "velocity-L2-error(rel)=" << error_abs/u_abs << endl;

  SuperL2Norm3D<T> diffL2NormObj(uSolLattice-u, superGeometry, *objectiveDomain);
  SuperL2Norm3D<T> absL2NormObj(uSolLattice, superGeometry, *objectiveDomain);
  diffL2NormObj(&error_abs, nullptr);
  absL2NormObj(&u_abs, nullptr);
  clout << "ObjInd: ";
  clout << "velocity-L2-error(abs)=" << error_abs << "; \t";
  clout << "velocity-L2-error(rel)=" << error_abs/u_abs << endl;

  // Plot L2-error
  static Gnuplot<T> gplot("L2error");
  gplot.setData(true, error_abs/u_abs, "Obj: velocity-L2-error(rel)");
  gplot.writePNG();
}


/// Solver
template<typename S, typename T, template<typename U> class Lattice>
void Solver3D<S,T,Lattice>::prepareGeometry()
{
  /*
   * - Material 0: Outside of this problem
   * - Material 1: inner flow field with porosity=1
   * - Material 2: Border cells
   * - Material 3: Poiseuille flow
   * - Material 4: Pressure boundary
   * - Material 5: Solid object from artificial data
   * - Material 6: Design area
   */
  clout << "Prepare Geometry ..." << std::endl;

  // Read mesh geometry from XML file
  IndicatorF3D<T>* cube = createIndicatorF3D<T>((*this->_xml)["Application"]["mesh"], false);
  // Create CuboidGeometry3D from mesh indicator
  _cGeometry = new CuboidGeometry3D<T>(*cube, this->_converter->getPhysDeltaX(), singleton::mpi().getSize() );
  // Create HeuristicLoadBalancer from CuboidGeometry3D
  _loadBalancer = new HeuristicLoadBalancer<T>(*_cGeometry);
  // Create SuperGeometry3D with __overlap
  _bg = new SuperGeometry3D<T>(*_cGeometry, *_loadBalancer, __overlap);

  objectiveDomain = new SuperIndicatorFfromIndicatorF3D<T>(*createIndicatorF3D<T>((*this->_xml)["ObjectiveDomain"], false), *_bg);

  // Material 2 everywhere
  _bg->rename(0,2, *cube);
  // Material 1 inside of border
  _bg->rename(2,1,1,1,1);

  ///adjusted to fit vti data
  Vector<T,3> center1,center2;
  center1 = { 66e-4, 54e-4, -this->_converter->getPhysDeltaX() };
  center2 = { 66e-4, 54e-4, +this->_converter->getPhysDeltaX() };

  IndicatorCylinder3D<T> inflow(center1,center2, 66e-4-this->_converter->getPhysDeltaX()*1.1);
  _bg->rename( 2,3,inflow );

  center1[1] += 296e-4; //length of domain
  center2[1] += 296e-4;
  IndicatorCylinder3D<T> outflow(center1,center2, 66e-4-this->_converter->getPhysDeltaX()*1.1);
  _bg->rename( 2,4,outflow );

  // Indicators for material 5 (minDesignDomain) and material 6 (designDomain)
  IndicatorF3D<T>* simulationObject = createIndicatorF3D<T>((*this->_xml)["SimulationObject"]);
  IndicatorF3D<T>* designDomain = createIndicatorF3D<T>((*this->_xml)["DesignDomain"]);

  if (optimizationStarted) {
    // Material 6: optimisation area
    _bg->rename(1, 6, *designDomain);
  } else {
    // Material 5: concrete object in simulation
    _bg->rename(1, 5, *simulationObject);
  }


//  _bg->clean(); --> dann funktioniert optimierung nicht mehr! warum??
  _bg->innerClean();
  _bg->checkForErrors();

  _bgStat = &_bg->getStatistics();
//  _bg->communicate(); --> nötig?

  _bg->getCuboidGeometry().print();
  _bg->getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}



template<typename S, typename T, template<typename U> class Lattice>
void Solver3D<S,T,Lattice>::prepareLattice()
{
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = this->_converter->getLatticeRelaxationFrequency();
  print(omega, "omega", clout);

  // material==0 -->do nothing
  _lattice->defineDynamics(*this->_bg, 0,  &instances::getNoDynamics<T, Lattice>());

  // material==1 -->bulk dynamics
  _lattice->defineDynamics(*this->_bg, 1, this->_bulkDynamics);

  // material==2 --> primal: bulk dynamics, dual: BounceBack
  if (!(this->_dualMode) ) {
    _lattice->defineDynamics(*this->_bg, 2, &instances::getBounceBack<T, Lattice>());
    for (int i: {3,4}) {
      _lattice->defineDynamics(*this->_bg, i, this->_bulkDynamics);
    }
  }
  if ((this->_dualMode) ) {
    for (int i: {2,3,4}) {
      _lattice->defineDynamics(*this->_bg, i, &instances::getBounceBack<T, Lattice>());
    }
  }

  if (optimizationStarted) {
    // material==6 --> bulk dynamics
    _lattice->defineDynamics(*this->_bg, 6, this->_bulkDynamics);
  } else {
    // material==5 --> bulkDynamics/BounceBack (BounceBack more realistic, bulkDyn easier to optimize against
    _lattice->defineDynamics(*this->_bg, 5, this->_bulkDynamics);
//    _lattice->defineDynamics(*this->_bg, 5, &instances::getBounceBack<T, Lattice>());
  }

  // Setting of the boundary conditions
  if (!(this->_dualMode) ) {
    _onBc->addVelocityBoundary(*this->_bg, 3, omega);
//    _onBc->addVelocityBoundary(*this->_bg, 4, omega);
    _onBc->addPressureBoundary(*this->_bg, 4, omega);
  }
  _lattice->initialize();

  clout << "Prepare Lattice... OK" << std::endl;
}


/// Set Boundary Values
template<typename S, typename T, template<typename U> class Lattice>
void Solver3D<S,T,Lattice>::setBoundaryValues(int iT)
{
  static OstreamManager clout("setBoundaryValues");

  int itStart = this->_converter->getLatticeTime(this->_startTime);

  if (iT==0) {

    AnalyticalConst3D<T,T> zero(0);
    AnalyticalConst3D<T,T> one(1);

    AnalyticalConst3D<T,T> rhoF(1);
    std::vector<T> velocity(3,T());
    AnalyticalConst3D<T,T> uF(velocity);

    for ( int i: {0,1,2,3,4} ) {
      _lattice->defineExternalField(*this->_bg, i, Lattice<T>::ExternalField::porosityIsAt,
                                    Lattice<T>::ExternalField::sizeOfPorosity, one);
    }

    if (optimizationStarted) {
      _lattice->defineExternalField(*this->_bg, 5, Lattice<T>::ExternalField::porosityIsAt,
                                    Lattice<T>::ExternalField::sizeOfPorosity, *(this->_startExternalField));
      _lattice->defineExternalField(*this->_bg, 6, Lattice<T>::ExternalField::porosityIsAt,
                                    Lattice<T>::ExternalField::sizeOfPorosity, *(this->_startExternalField));
    } else {
      _lattice->defineExternalField(*this->_bg, 5, Lattice<T>::ExternalField::porosityIsAt,
                                    Lattice<T>::ExternalField::sizeOfPorosity, zero);
    }

    for ( int i: {1,2,3,4,5,6} ) {
      _lattice->iniEquilibrium(*(this->_bg), i, rhoF, uF);
    }
  }

  if (iT <= itStart) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( itStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity;

    if (iT%50==0  && !(this->_dualMode)) {
      maxVelocity = 2.25*frac[0]*this->_converter->getCharLatticeVelocity();
      CirclePoiseuille3D<T> poiseuilleU(*this->_bg, 3, maxVelocity);
      _lattice->defineU( *this->_bg, 3, poiseuilleU );
    }
  }
}


/// Opti: Objective functions
template<typename T, template<typename U> class Lattice>
bool Objective3D<T,Lattice>::operator()(T output[], const int x[])
{
  static OstreamManager clout("Objective3D");

  error(this->_geometry, this->_sLattice, this->_converter );

  SuperLatticePhysVelocity3D<T,Lattice> physLatticeVelocity(this->_sLattice,this->_converter);
  RelativeDifferenceObjective3D<T,Lattice> objective(this->_sLattice, physLatticeVelocity, *solution, this->_geometry, *objectiveDomain);

  return objective(output,x);
}


template<typename T, template<typename U> class Lattice>
DObjectiveDf3D<T,Lattice>::DObjectiveDf3D(SuperLattice3D<T,Lattice>& sLattice, SuperGeometry3D<T>& geometry,
    Controller<T>& controller, UnitConverter<T,Lattice>& converter)
  : ControlledSuperLatticeF3D<T,Lattice>(sLattice, geometry, controller, converter, Lattice<T>::q)
{
  this->getName() = "dObjectivedF";

  static SuperLatticePhysVelocity3D<T,Lattice> physLatticeVelocity(this->_sLattice,this->_converter);
  static SuperLatticeDphysVelocityDf3D<T,Lattice> dPhysVelocityf(this->_sLattice, this->_converter);
  
  _dObjectiveDf = new DrelativeDifferenceObjectiveDf3D<T,Lattice>(this->_sLattice, physLatticeVelocity, dPhysVelocityf, *solution,
      this->_geometry, *objectiveDomain);
}


template<typename T, template<typename U> class Lattice>
bool DObjectiveDf3D<T,Lattice>::operator()(T output[], const int latticeR[])
{
  static OstreamManager clout(std::cout,"DObjectiveDf3D");

  return (*_dObjectiveDf)(output,latticeR);
}


template<typename T, template<typename U> class Lattice>
bool DObjectiveDControl3D<T,Lattice>::operator()(T output[], const int x[])
{
  static OstreamManager clout(std::cout,"DObjectiveDControl3D");
  // Objective does not depend on control! -> deriv is 0
  output[0] = T();
  output[1] = T();
  output[2] = T();
  return true;
}



/// Get results
template<typename S, typename T, template<typename U> class Lattice>
void Solver3D<S,T,Lattice>::getResults(int iT)
{
  // Open Output Files
  SuperVTMwriter3D<T> simuWriter("Simulation_FlowField");
  SuperVTMwriter3D<T> optiWriter("Optimization_Porosity");
  SuperVTMwriter3D<T> errorWriter("Error_FlowField");

  // Optimization Mode
  if ((this->_dualMode) ) {
    if (iT==0) {
      if (counterOptiStep==0) {
        optiWriter.createMasterFile();
      }

      // ---------- add porosity ---------- //
      SuperExternal3D<T, Lattice> exField(*(this->_bg), *_lattice, Lattice<T>::ExternalField::porosityIsAt,
                                          Lattice<T>::ExternalField::sizeOfPorosity, __overlap);
      exField.communicate();

      SuperLatticePorosity3D<T, Lattice> d(*_lattice);

      // Write Porosity into VTK File
      optiWriter.addFunctor( d );
      optiWriter.write(counterOptiStep);

      SuperLatticeGeometry3D<T, Lattice> geometry(*_lattice, *(this->_bg));
      optiWriter.write(geometry,1);

      counterOptiStep++;
    }
  }

  // Simulation Mode
  if (!(this->_dualMode) ) {
    SuperVTMwriter3D<T>* writer = optimizationStarted ? &errorWriter : &simuWriter;
    if ( iT == 0 ) {
      if (!(optimizationStarted)) {
        // Writes the geometry, cuboid no. and rank no. as vti file for visualization
        SuperLatticeGeometry3D<T, Lattice> geometry(*_lattice, *(this->_bg));
        SuperLatticeCuboid3D<T, Lattice> cuboid(*_lattice);
        SuperLatticeRank3D<T, Lattice> rank(*_lattice);

        writer->write(cuboid);
        writer->write(geometry);
        writer->write(rank);
        writer->createMasterFile();
      } else {
        if (counterOptiStep == 0) {
          writer->createMasterFile();
        } else {
          numberOfSimulations[counterOptiStep]++;
          if (singleton::mpi().isMainProcessor()) {
            std::ofstream plotOut(errorPlotFile, std::ios::trunc);
            for (unsigned int i=0; i<errorNorm.size(); i++) {
              plotOut << std::setprecision(16) << i << "\t" << errorNorm[i] << "\t" << errorNormInside[i] << "\t" << numberOfSimulations[i] << endl;
            }
            plotOut.close();
          }
        }

        // add new plot entries
        if (errorNorm.size() < counterOptiStep + 1) {
          errorNorm.push_back(0);
          errorNormInside.push_back(0);
          numberOfSimulations.push_back(1);
        }
      }
    }

    // Opti mode: Overwrite flow field for iT again and again
    if (optimizationStarted) {
      // get error functor
      SuperLatticePhysVelocity3D<T, Lattice> velocity(*_lattice, *(this->_converter));
      SuperLatticeFfromAnalyticalF3D<T, Lattice> solVelocity_orig(*solution, *_lattice);
      SuperExtractIndicatorF3D<T> solVelocity(solVelocity_orig, *objectiveDomain);

      // Write Velocity and Pressure into Simulation File
      SuperLatticePhysVelocity3D<T, Lattice> sVelocity(*_lattice, *(this->_converter));
      SuperLatticeFfromAnalyticalF3D<T,Lattice> mriField_orig(*solution, *_lattice);
      SuperExtractIndicatorF3D<T> mriField(mriField_orig, *objectiveDomain);

      writer->addFunctor(sVelocity);
      writer->addFunctor(mriField);
      writer->write(counterOptiStep);

      // ---------- plot L2 error ------------//
      T error_abs, vel_abs, error_abs_inside, vel_abs_inside;

      SuperL2Norm3D<T> diffL2Norm(solVelocity - velocity, *_bg, *objectiveDomain);
      SuperL2Norm3D<T> diffL2NormInside(solVelocity - velocity, *_bg, 6);
      SuperL2Norm3D<T> absL2Norm(solVelocity, *_bg, *objectiveDomain);
      SuperL2Norm3D<T> absL2NormInside(solVelocity, *_bg, 6);

      diffL2Norm(&error_abs, nullptr);
      diffL2NormInside(&error_abs_inside, nullptr);
      absL2Norm(&vel_abs, nullptr);
      absL2NormInside(&vel_abs_inside, nullptr);

      errorNorm[counterOptiStep] = error_abs / vel_abs;
      errorNormInside[counterOptiStep] = error_abs_inside / vel_abs_inside;

      errorNormCounter++;
    } else {
      // Write Velocity and Pressure into Simulation File
      SuperLatticePhysVelocity3D<T, Lattice> velocity(*_lattice, *(this->_converter));
      SuperLatticePhysPressure3D<T, Lattice> pressure(*_lattice, *(this->_converter));

      writer->addFunctor(velocity);
      writer->addFunctor(pressure);
      writer->write(iT);
    }
  }
}




/*********************************************/
/*                 MAIN                      */
/*********************************************/
int main(int argc, char *argv[])
{
  // initialize OpenLB and some of its modules
  olbInit(&argc, &argv);
  OstreamManager clout(std::cout,"Main");

  std::string opti_xml = "parameter.xml";
  clout << "Reading XML File \"" << opti_xml << "\"..." << endl;
  XMLreader opti_config(opti_xml);

  clout.setMultiOutput(opti_config["Output"]["MultiOutput"].get<bool>());
  output_dir = "output/";/*
               + opti_config["Application"]["OutputName"].get<string>() + "_"
               + opti_config["Optimization"]["StartValue"].get<string>() + "_"
               + opti_config["Application"]["Discretization"]["Resolution"].get<string>() + "/";*/
  singleton::directories().setOutputDir(output_dir);

  /*UnitConverter<S,Lattice>* lbc;
  SuperLattice3D<S, Lattice>* sl;

  /// === Simulate Artificial Data ===
  XMLreader solver_config(opti_xml);

  // initialize the solver and the problem from xml
  Solver3D<S, S, Lattice> solverProblem(&solver_config);

  clout << endl << "=== Simulate Data ===" << endl << endl;
  // init solver3D
  solverProblem.init();
  singleton::directories().setOutputDir(output_dir);

  // solve primal problem to get u* and p* (reference vel.- and press.-field)
  solverProblem.solve();

  clout << "Simulation finished" << endl << endl;

  // Fill objects pointers
  lbc = solverProblem.getConverter();
  lbc->print();
  sl = solverProblem.getLattice();
  clout << "Creating SuperLatticePhysVelocity3D..." << endl;
  SuperLatticePhysVelocity3D<S, Lattice> velocityLattice(*sl,*lbc);

  // Store the Analytical functor of the velocity functor in the global variable solution
  clout << "Interpolating over velocity functor..." << endl;
  solution = new AnalyticalFfromSuperF3D<T>(velocityLattice, false, __overlap, false);

  clout << "Interpolation finished. Starting optimisation on Analytical Functor..." << endl << endl;

  SuperVTMwriter3D<T> solutionWriter("solution");
  solutionWriter.addFunctor(velocityLattice);
  solutionWriter.write();
*/
  
  // === Read MRI Data === //
  std::string dataName = "result.vti"; //opti_config["Application"]["MRIdata"].get<string>()+".vti";

  BlockVTIreader3D<T,BaseType> readerData (dataName, "physVelocity");
  BlockData3D<T,BaseType> data = readerData.getBlockData();
  BlockDataF3D<T, BaseType> velocityF(data);
  BlockExtractIndicatorF3D<T> velocityLattice(velocityF, *objectiveDomain);

  // Store the Analytical Functor of the velocity functor in the global variable solution
  clout << "Interpolating over velocity functor..." << endl;
  // Read Spacing_
  XMLreader MRI_config(dataName);

//velocityLattice,readerData...
  solution = new SpecialAnalyticalFfromBlockF3D<T,BaseType>(velocityLattice, readerData.getCuboid(), Vector<T,3> (2e-4, 2e-4, 2e-4));

  clout << "Interpolation finished." << endl << endl;
  

  // === Solve Optimization Problem ===
  clout << "=== Start Optimization ===" << endl << endl;
  errorPlotFile = "errorPlot_.dat";/*
                  + opti_config["Application"]["OutputName"].get<string>() + "_"
                  + opti_config["Optimization"]["StartValue"].get<string>() + "_"
                  + opti_config["Application"]["Discretization"]["Resolution"].get<string>() + ".dat";*/

  // Set optimizationStarted variable (for output)
  optimizationStarted = true;


  // Solve optimization problem
  clout << "== Initialize problem ==" << std::endl;
  OptiSolver<S, T, Lattice> optiProblem(&opti_config);
  singleton::directories().setOutputDir(output_dir);
  clout << "== Solve Problem ==" << std::endl;
  optiProblem.solve();
  clout << "Optimization finished." << std::endl;

  return 0;
}
