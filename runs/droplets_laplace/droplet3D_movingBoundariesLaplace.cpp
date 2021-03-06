/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

/* rayleighTaylor3d.cpp:
 * Rayleigh-Taylor instability in 3D, generated by a heavy
 * fluid penetrating a light one. The multi-component fluid model
 * by X. Shan and H. Chen is used. This example shows the usage
 * of multicomponent flow and periodic boundaries.
 */


#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!
#include <cstdlib>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR ShanChenDynOmegaForcedD3Q19Descriptor


// Parameters for the simulation setup
const double nx   = 0.4; // in m
const double ny   = 0.1;
const double nz   = 0.1;
const double radius = 12/10. * 0.02;
const int maxIter  = 3000;

const int interval = 20; // in s
const T epsilon = 1e-5;

const int N = 160;

const double G = 8.;

T movingWallVelocity = 0.0;


// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry3D<T>& superGeometry, UnitConverter<T, DESCRIPTOR> const& converter) {

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename(0,1);

  std::vector<T> origin4(3,T());
    origin4[0] = 0;//-2.*converter.getLatticeL();
    origin4[1] = 0;//-2.*converter.getLatticeL();
    origin4[2] = 0;//-2.*converter.getLatticeL();


    std::vector<T> origin3(3,T());
    origin3[0] = 0;
    origin3[1] = ny;
    origin3[2] = 0;

    Vector<T,3> origin5( 0,0,0 );
    Vector<T,3> origin6( 0,0,nz);



    std::vector<T> extend4(3,T());
    extend4[0] = nx;
    extend4[1] = 0;
    extend4[2] = nz;

    std::vector<T> extend3(3,T());
    extend3[0] = nx;
    extend3[1] = 0;
    extend3[2] = nz;








    IndicatorCuboid3D<T> bottom(extend4, origin4);
    IndicatorCuboid3D<T> top(extend3, origin3);
    //IndicatorCuboid3D<T> front(extend5, origin5);
    //IndicatorCuboid3D<T> back(extend6, origin6);



    superGeometry.rename(1,3,top);
    superGeometry.rename(1,4,bottom);
    //superGeometry.rename(1,5,front);
    //superGeometry.rename(1,6,back);

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;

}


void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& sLatticeOne,
                     SuperLattice3D<T, DESCRIPTOR>& sLatticeTwo,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics1,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics2,
                     Dynamics<T, DESCRIPTOR>& bounceBackRho0,
                     Dynamics<T, DESCRIPTOR>& bounceBackRho1,
                     SuperGeometry3D<T>& superGeometry,
                     sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> &bcOne,
                     sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> &bcTwo,
                     UnitConverter<T, DESCRIPTOR> const& converter) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // define lattice Dynamics
  sLatticeOne.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLatticeTwo.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLatticeOne.defineDynamics( superGeometry, 1, &bulkDynamics1 );
  sLatticeTwo.defineDynamics( superGeometry, 1, &bulkDynamics2 );

  sLatticeOne.defineDynamics( superGeometry, 3, &bulkDynamics1 );
  sLatticeTwo.defineDynamics( superGeometry, 3, &bulkDynamics2 );



  sLatticeOne.defineDynamics( superGeometry, 4, &bulkDynamics1 );
  sLatticeTwo.defineDynamics( superGeometry, 4, &bulkDynamics2 );

  bcOne.addVelocityBoundary( superGeometry, 3, converter.getLatticeRelaxationFrequency());
  bcOne.addVelocityBoundary( superGeometry, 4, converter.getLatticeRelaxationFrequency());


  bcTwo.addVelocityBoundary( superGeometry, 3, converter.getLatticeRelaxationFrequency());
  bcTwo.addVelocityBoundary( superGeometry, 4, converter.getLatticeRelaxationFrequency());


  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLatticeOne,
                        SuperLattice3D<T, DESCRIPTOR>& sLatticeTwo,
                        T force, int iT, SuperGeometry3D<T>& superGeometry ) {

  if ( iT==0 ) {

      std::vector<T> v(3,T());
      AnalyticalConst3D<T,T> zeroV(v);
      AnalyticalConst3D<T,T> zero(1e-4);

      AnalyticalConst3D<T,T> rhoOne(1.0);
      AnalyticalConst3D<T,T> rhoNearZero(1e-5);

      std::vector<T> center = {nx/2, ny/2, nz/2};
      std::vector<T> center1 = {nx/2-radius, ny/2,nz/2};
      std::vector<T> center2 = {nx/2+radius, ny/2,nz/2};

      std::cout << center1[0] << " " << center1[1] << " " << center1[2] << std::endl;
      std::cout << center2[0] << " " << center2[1] << " " << center2[2] << std::endl;

      //SmoothIndicatorCylinder3D<T,T> sphere(center1,center2,radius, radius*0.1);
      SmoothIndicatorSphere3D<T,T> sphere(center, radius, radius*0.2, 1);

      AnalyticalConst3D<T,T> one(1.0);

      AnalyticalIdentity3D<T,T> rhoLatticeOne((rhoOne - rhoNearZero) * sphere + rhoNearZero);
      AnalyticalIdentity3D<T,T> rhoLatticeTwo((rhoOne - rhoNearZero) * (one - sphere) + rhoNearZero);


      std::vector<T> F(3,T());
      F[1] = -force;
      AnalyticalConst3D<T,T> f(F);

    // for each material set the defineRhou and the Equilibrium

    sLatticeOne.defineRhoU( superGeometry, 1, rhoLatticeOne, zeroV );
    sLatticeOne.iniEquilibrium( superGeometry, 1, rhoLatticeOne, zeroV );
    sLatticeOne.defineExternalField( superGeometry, 1,
                                     DESCRIPTOR<T>::ExternalField::externalForceBeginsAt,
                                     DESCRIPTOR<T>::ExternalField::sizeOfExternalForce, f );
    sLatticeTwo.defineRhoU( superGeometry, 1, rhoLatticeTwo, zeroV );
    sLatticeTwo.iniEquilibrium( superGeometry, 1, rhoLatticeTwo, zeroV );
    sLatticeTwo.defineExternalField( superGeometry, 1,
                                     DESCRIPTOR<T>::ExternalField::externalForceBeginsAt,
                                     DESCRIPTOR<T>::ExternalField::sizeOfExternalForce, f );

    //T movingWallVelocity = 0.;

    std::vector<T> u = {movingWallVelocity,0,0};
    AnalyticalConst3D<T,T> vel(u);
    std::vector<T> u2 = {-movingWallVelocity,0,0};
    AnalyticalConst3D<T,T> vel2(u2);

    sLatticeTwo.defineRhoU( superGeometry, 3, rhoLatticeTwo, vel );
    sLatticeTwo.iniEquilibrium( superGeometry, 3, rhoLatticeTwo, vel );
    sLatticeTwo.defineRhoU( superGeometry, 4, rhoLatticeTwo, vel2 );
    sLatticeTwo.iniEquilibrium( superGeometry, 4, rhoLatticeTwo, vel2 );

    sLatticeOne.defineRhoU( superGeometry, 3, rhoLatticeOne, vel );
    sLatticeOne.iniEquilibrium( superGeometry, 3, rhoLatticeOne, vel );
    sLatticeOne.defineRhoU( superGeometry, 4, rhoLatticeOne, vel2 );
    sLatticeOne.iniEquilibrium( superGeometry, 4, rhoLatticeOne, vel2 );





    // Make the lattice ready for simulation
    sLatticeOne.initialize();
    sLatticeTwo.initialize();
  }
}

void getResults( SuperLattice3D<T, DESCRIPTOR>& sLatticeTwo,
                 SuperLattice3D<T, DESCRIPTOR>& sLatticeOne, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> const& converter) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "droplet3D" );

  const int vtkIter  = converter.getLatticeTime(15.);

  const int statIter = converter.getLatticeTime(15.);

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLatticeOne, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLatticeOne );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLatticeOne );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

  
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter==0 && iT > 0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    sLatticeOne.getStatistics().print(iT,iT);

    clout << "averageRhoFluidOne="   << sLatticeOne.getStatistics().getAverageRho();
    clout << "; averageRhoFluidTwo=" << sLatticeTwo.getStatistics().getAverageRho() << std::endl;
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    clout << "Writing VTK ..." << std::endl;
    SuperLatticeVelocity3D<T, DESCRIPTOR> velocity( sLatticeOne );
    SuperLatticeDensity3D<T, DESCRIPTOR> density( sLatticeOne );
    SuperLatticeGeometry3D<T, DESCRIPTOR> materialOne(sLatticeOne,superGeometry);
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( density );
    vtmWriter.addFunctor( materialOne );



    SuperLatticeVelocity3D<T, DESCRIPTOR> velocity2( sLatticeTwo );
    velocity2.getName() = "velocity2";
    SuperLatticeDensity3D<T, DESCRIPTOR> density2( sLatticeTwo );
    density2.getName() = "density2";    

    vtmWriter.addFunctor( velocity2 );
    vtmWriter.addFunctor( density2 );
    
    vtmWriter.write( iT );

    BlockReduction3D2D<T> planeReduction( density, {0, 0, -1});
    BlockGifWriter<T> gifWriter;
    gifWriter.write( planeReduction, iT, "density" );

    clout << "Writing VTK ... OK" << std::endl;
  }


  // Saves lattice data
  if ( iT == maxIter -1 ) {
    //clout << "Checkpointing the system at t=" << iT << endl;
    //sLatticeOne.save( "latticeOne.checkpoint" );
    //sLatticeTwo.save( "latticeTwo.checkpoint" );
  }
  /*
  if ( iT == 0 ) {
    clout << "Checkpointing the system at t=" << iT << endl;
    sLatticeOne.load( "latticeOne.checkpoint" );
    sLatticeTwo.load( "latticeTwo.checkpoint" );
  }*/

}


int main( int argc, char *argv[] ) {
  T name = radius/0.0025;
  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp_G8_r_"+std::to_string(name)+"/");//movingWallVelocity)+"/" );
  OstreamManager clout( std::cout,"main" );

  T force        = 0./(T)ny/(T)ny;


  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      int {nx/(double)N},                        // resolution: number of voxels per charPhysL
      (T)   0.98,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)   nx,       // charPhysLength: reference length of simulation geometry
      (T)   0.001,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)   1e-6, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)   998.0                       // physDensity: physical density in __kg / m^3__
    );


  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights

  std::vector<T> origin = {0., 0., 0.};
  std::vector<T> extend = {nx, ny, nz};
  IndicatorCuboid3D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cGeometry( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidGeometry3D<T> cGeometry( cuboid, converter.getPhysDeltaX(), 1 );
#endif

  cGeometry.setPeriodicity( true, false, true );

  HeuristicLoadBalancer<T> loadBalancer( cGeometry );

  SuperGeometry3D<T> superGeometry( cGeometry,loadBalancer,2 );

  prepareGeometry( superGeometry , converter);

  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLatticeOne( superGeometry );
  SuperLattice3D<T, DESCRIPTOR> sLatticeTwo( superGeometry );

  ForcedBGKdynamics<T, DESCRIPTOR> bulkDynamics1 (
    converter.getLatticeRelaxationFrequency(), instances::getExternalVelocityMomenta<T,DESCRIPTOR>() );
  ForcedBGKdynamics<T, DESCRIPTOR> bulkDynamics2 (
    converter.getLatticeRelaxationFrequency(), instances::getExternalVelocityMomenta<T,DESCRIPTOR>() );


  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondOne(sLatticeOne);
  createLocalBoundaryCondition3D<T, DESCRIPTOR>(sBoundaryCondOne);

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondTwo(sLatticeTwo);
   createLocalBoundaryCondition3D<T, DESCRIPTOR>(sBoundaryCondTwo);

  // A bounce-back node with fictitious density 1,
  //   which is experienced by the partner fluid
  BounceBack<T, DESCRIPTOR> bounceBackRho1( 1. );
  // A bounce-back node with fictitious density 0,
  //   which is experienced by the partner fluid
  BounceBack<T, DESCRIPTOR> bounceBackRho0( 0. );

  std::vector<T> rho0;
  rho0.push_back( 1 );
  rho0.push_back( 1 );
  ShanChen93<T,T> interactionPotential;
  ShanChenForcedGenerator3D<T,DESCRIPTOR> coupling( G,rho0,interactionPotential );

  sLatticeOne.addLatticeCoupling( superGeometry, 1, coupling, sLatticeTwo );

  prepareLattice( sLatticeOne, sLatticeTwo, bulkDynamics1, bulkDynamics2,
                  bounceBackRho0, bounceBackRho1, superGeometry, sBoundaryCondOne, sBoundaryCondTwo, converter );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << endl;
  Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge( converter.getLatticeTime(interval), epsilon );

  for ( iT=0; iT<maxIter; ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLatticeOne, sLatticeTwo, force, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLatticeOne.collideAndStream();
    sLatticeTwo.collideAndStream();

    sLatticeOne.communicate();
    sLatticeTwo.communicate();

    sLatticeOne.executeCoupling();
    //sLatticeTwo.executeCoupling();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLatticeTwo, sLatticeOne, iT, superGeometry, timer, converter );

    /*clout << "averageRhoOne= " << sLatticeOne.getStatistics().getAverageRho();
    clout << ";  averageRhoTwo= " << sLatticeTwo.getStatistics().getAverageRho() << endl;
    clout << "maxU one= " << sLatticeOne.getStatistics().getMaxU();
    clout << ";  maxU Two= " << sLatticeTwo.getStatistics().getMaxU() << std::endl;

*/

	if (sLatticeOne.getStatistics().getAverageRho() != sLatticeOne.getStatistics().getAverageRho())
	{
		std::cout << "\n Rho diverges ... simulation killed \n" << std::endl;
		break;


	}






    if(converge.hasConverged()) break;

    converge.takeValue( sLatticeOne.getStatistics().getAverageEnergy(), true );

  }
  std::cout << "converged or ended normally =) " << std::endl;
  timer.stop();
  timer.printSummary();
}

