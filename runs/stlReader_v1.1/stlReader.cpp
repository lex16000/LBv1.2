/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
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

/* cylinder3d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 * It also shows the usage of the STL-reader and explains how
 * to set boundary conditions automatically.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor


// Parameters for the simulation setup
const int N = 1;        // resolution of the model
const int M = 1;        // time discretization refinement
     // Reynolds number
const T maxPhysT = 10.; // max. simulation time in s, SI unit


T radiusInlet =  0.50; // m
T radiusOutlet = 0.50; // m

const int interval = 2; // in s
const T epsilon = 1e-4;


// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( LBconverter<T> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,1,1,1);

  //superGeometry.rename( 2,1,stlReader );
  superGeometry.clean();

  Vector<T,3> origin = superGeometry.getStatistics().getMinPhysR( 2 );




  origin[0] += converter.getLatticeL()/2.;
  origin[2] += converter.getLatticeL()/2.;

  Vector<T,3> extend = superGeometry.getStatistics().getMaxPhysR( 2 );
  extend[0] = extend[0]-origin[0]-converter.getLatticeL()/2.;
  extend[2] = extend[2]-origin[2]-converter.getLatticeL()/2.;

  

  // Set material number for inflow
  origin[1] = superGeometry.getStatistics().getMinPhysR( 2 )[1]-converter.getLatticeL();
  extend[1] = 2*converter.getLatticeL();
  IndicatorCuboid3D<T> inflow( extend,origin );
  superGeometry.rename( 2,3,inflow );



  // Set material number for outflow
  origin[1] = superGeometry.getStatistics().getMaxPhysR( 2 )[1]-converter.getLatticeL();
  extend[1] = 2*converter.getLatticeL();
  IndicatorCuboid3D<T> outflow( extend,origin );
  superGeometry.rename( 2,4,outflow );




  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();
  clout << "; press enter" << std::endl;
  std::cin.ignore(10, '\n');
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     LBconverter<T> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                     STLreader<T>& stlReader,
                     SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getOmega();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  // Material=3 -->bulk dynamics (inflow)
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );

  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );
  //sLattice.defineDynamics( superGeometry, 5, &bulkDynamics );

  // Setting of the boundary conditions

  bc.addVelocityBoundary( superGeometry, 3, omega );
  bc.addPressureBoundary( superGeometry, 4, omega );




  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( superGeometry, 1, rhoF, uF );
  sLattice.iniEquilibrium( superGeometry, 1, rhoF, uF );
  sLattice.defineRhoU( superGeometry, 3, rhoF, uF );
  sLattice.iniEquilibrium( superGeometry, 3, rhoF, uF );
  sLattice.defineRhoU( superGeometry, 4, rhoF, uF );
  sLattice.iniEquilibrium( superGeometry, 4, rhoF, uF );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        LBconverter<T> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.numTimeSteps( maxPhysT*0.1 );
  int iTupdate = 30;

  if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[1] = 2.25*frac[0]*converter.getLatticeU();

    T distance2Wall = converter.getLatticeL()/2.;
    //RectanglePoiseuille3D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    AnalyticalConst3D<T,T> uInlet(0,frac[0]*converter.getLatticeU(),0);
    //sLattice.defineU( superGeometry, 3, poiseuilleU );
    sLattice.defineU( superGeometry, 3, uInlet);

    clout << "step=" << iT << "; maxVel=" << maxVelocity[1] << std::endl;
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 LBconverter<T>& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer,
                 STLreader<T>& stlReader ) {

  OstreamManager clout( std::cout,"getResults" );

  T lengthOfPipe = superGeometry.getStatistics().getMaxPhysR( 4 )[1];

  SuperVTMwriter3D<T> vtmWriter( "cylinder3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeDensity3D<T, DESCRIPTOR> density( sLattice );
  SuperLatticeGeometry3D<T, DESCRIPTOR> material( sLattice, superGeometry );
  //SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, superGeometry, stlReader, 5 );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( density );
  vtmWriter.addFunctor( material );
  //vtmWriter.addFunctor( yPlus );

  const int vtkIter  = 30;//converter.numTimeSteps( maxPhysT/25 );
  const int statIter = 30;//converter.numTimeSteps( maxPhysT/25 );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%vtkIter == 0 ) {
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockLatticeReduction3D<T, DESCRIPTOR> planeReduction( normVel, 1, 0, 0 );
    BlockGifWriter<T> gifWriter;
    //gifWriter.write(planeReduction, 0, 0.7, iT, "vel"); //static scale
    gifWriter.write( planeReduction, iT, "vel" ); // scaled

    cout << "writing vtk... " << (double)iT/converter.numTimeSteps(maxPhysT)*100 << " %" << endl;
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.physTime( iT ) );
 

   

  }
}

int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);
T Re = 1e-3;
  LBconverter<T> converter(
	     3,                                     // dim
	     ( T ) 0.5,                     // latticeL_
	     ( T ) 0.0001,                          // latticeU_
         ( T ) 2e-5,                     // charNu_
	     ( T ) 90,                         // charL_ = 1,
         ( T ) 0.01,                            // charU_ = 1,
         ( T ) 1.2,                            // charRho_ = 1,
	     ( T ) 0                                // charPressure_ = 0
	   );
	   converter.print();
	   writeLogFile( converter, "porousBlock" );




  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  STLreader<T> stlReader( "USS_Constitution.stl", converter.getLatticeL(), 1 );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getLatticeL() );

  clout << stlReader.getMin()[0] << " " << stlReader.getMin()[1] << " " << stlReader.getMin()[2] << std::endl;
  clout << stlReader.getMax()[0] << " " << stlReader.getMax()[1] << " " << stlReader.getMax()[2] << std::endl;


  clout << "no timesteps: " << converter.numTimeSteps(maxPhysT) << std::endl;
  std::cout << " enter " << std::endl;
  cin.ignore(10, '\n');


  std::vector<T> origin = {-16.6, -45.3, -0.155};
  std::vector<T> extend = {28.5, 97.4,66.1};
  IndicatorCuboid3D<T> cuboidBig(extend, origin);

  IndicatorIdentity3D<T> domainIndicator( cuboidBig - stlReader );
 

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry3D<T> cuboidGeometry( domainIndicator, converter.getLatticeL(), noOfCuboids );
  //CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getLatticeL(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, domainIndicator, stlReader, superGeometry );
  //prepareGeometry( converter, extendedDomain, stlReader, superGeometry );



  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );
  BGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getOmega(), instances::getBulkMomenta<T, DESCRIPTOR>() );

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  //createInterpBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryCondition );
  createLocalBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );

  prepareLattice( sLattice, converter, bulkDynamics, sBoundaryCondition, sOffBoundaryCondition, stlReader, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer( converter.numTimeSteps( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge( converter.numTimeSteps(interval), epsilon );
clout << "no timesteps: " << converter.numTimeSteps(maxPhysT) << std::endl;
  for ( int iT = 0; iT < converter.numTimeSteps( maxPhysT ); ++iT ) {

if (iT == 0)
		{
			//std::cout << "loading checkpoint ..." << std::endl;
			//sLattice.load("sLattice.checkpointAb300000");
			//std::cout << "loading checkpoint ... OK" << std::endl;
		}
clout << "iT = " << iT << std::endl;
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, stlReader );

	if(converge.hasConverged()) {
	std::cout << "converged" << std::endl;	
	break;
	}

    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
    if(iT == converter.numTimeSteps( maxPhysT ) -1) clout << "converged or ended normally" << endl;
  }

  timer.stop();
  timer.printSummary();
}

