#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code;
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;
#define DESCRIPTOR D3Q19Descriptor


//simulation parameters
const bool bouzidiOn = true; // choice of boundary condition
const T maxPhysT = 15.;


const int nx = 62;
const int ny = 60;
const int nz = 148;

const T deltaX = 2e-4; // in m
const T lx = nx * deltaX;
const T ly = ny * deltaX;
const T lz = nz * deltaX;



const int interval = 1; // in s
const T epsilon = 4e-5;



void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& lattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& bc,
                     sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                     SuperGeometry3D<T>& superGeometry,
					 UnitConverter<T, DESCRIPTOR> const& converter) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // material=0 --> do nothing
  lattice.defineDynamics( superGeometry,0,&instances::getNoDynamics<T, DESCRIPTOR>() );

  // material=1 --> bulk dynamics
  lattice.defineDynamics( superGeometry,1,&bulkDynamics );


  // material=2 --> bounceBack dynamics
   lattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );




    // material=3 --> bulk dynamics + velocity (inflow)
    lattice.defineDynamics( superGeometry,3,&bulkDynamics );
    bc.addVelocityBoundary( superGeometry,3,omega );


  // material=4,5 --> bulk dynamics + pressure (outflow)
  lattice.defineDynamics( superGeometry,4,&bulkDynamics );
  bc.addPressureBoundary( superGeometry,4,omega );


  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  lattice.defineRhoU( superGeometry,1,rhoF,uF );
  lattice.iniEquilibrium( superGeometry,1,rhoF,uF );
  lattice.defineRhoU( superGeometry,3,rhoF,uF );
  lattice.iniEquilibrium( superGeometry,3,rhoF,uF );
  lattice.defineRhoU( superGeometry,4,rhoF,uF );
  lattice.iniEquilibrium( superGeometry,4,rhoF,uF );


  // Lattice initialize
  lattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}





void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
		UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.2 );
  int iTupdate = 30;

  if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    AnalyticalConst3D<T,T> uInlet(0,0,frac[0]*converter.getLatticeVelocity(converter.getCharPhysVelocity() ));
    //sLattice.defineU( superGeometry, 3, uInlet);
    //clout << "step=" << iT << "; maxVel=" << frac[0]*converter.getLatticeU() << std::endl;


    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[2] = 2.25*frac[0]*converter.getLatticeVelocity(converter.getCharPhysVelocity() );

    T distance2Wall = converter.getPhysDeltaX()/2.;
    RectanglePoiseuille3D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    sLattice.defineU( superGeometry, 3, poiseuilleU );
  }
}



void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 int iT,
                 SuperGeometry3D<T>& superGeometry,
                 UnitConverter<T, DESCRIPTOR> const& converter) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "read_voxel" );


  if ( iT==0 ) {
       SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
       SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
       SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
       vtmWriter.write( geometry );
       vtmWriter.write( cuboid );
       vtmWriter.write( rank );
       vtmWriter.createMasterFile();
     }


  SuperLatticeVelocity3D<T, DESCRIPTOR> latticeVelocity( sLattice );
  SuperLatticeDensity3D<T, DESCRIPTOR> latticeDensity( sLattice );
  SuperLatticeGeometry3D<T, DESCRIPTOR> material( sLattice, superGeometry );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> physVelocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeStrainRate3D<T, DESCRIPTOR> strain(sLattice,converter);
  vtmWriter.addFunctor( latticeVelocity );
  vtmWriter.addFunctor( latticeDensity );
  vtmWriter.addFunctor( material );
  vtmWriter.addFunctor( physVelocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( strain );

  const int statIter = converter.getLatticeTime( maxPhysT )/50.;
  const int checkpoint = statIter;



  // Writes output on the console
  if ( iT%statIter == 0 )
  {
	vtmWriter.write( iT );
	sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( latticeVelocity );
    BlockReduction3D2D<T> planeReduction( normVel, {0, 1, 0});
    heatmap::write(planeReduction, iT);


    //BlockGifWriter<T> gifWriter;
    //gifWriter.write( planeReduction, iT, "vel" ); // scaled

    cout << "writing vtk... " << (double)iT/converter.getLatticeTime(maxPhysT)*100 << " %" << endl;
  


  }


   if (iT%checkpoint == 0)
   {
    //clout << "loading the system" << endl;
    //sLatticeOne.load( "latticeOne.checkpoint" );
    //clout << "Checkpointing the system at timestep=" << iT << endl;
    //sLattice.save( "sLattice.checkpoint" );
   }



}


int main( int argc, char* argv[] ) {



  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );


  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      int {nx},                        // resolution: number of voxels per charPhysL
      (T)   0.56,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)   nx*deltaX,       // charPhysLength: reference length of simulation geometry
      (T)   0.001,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)   1e-6, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)   998.0                       // physDensity: physical density in __kg / m^3__
    );
//
//  LBconverter<T> converter(
//    ( int ) 3,                             // dim
//    ( T )   deltaX,                        // latticeL_
//    ( T )   0.005,                        // latticeU_
//    ( T )   1e-6,                // charNu_
//    ( T )   0.009,                           // charL_ = 1
//    ( T )   1e-3,                            // charU_ = 1
//    ( T )   998 // charRho = 998
//  );
  converter.print();
  converter.write("simpleTube" );


  // === 2nd Step: Prepare Geometry ===



  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  //const int noOfCuboids = 1;

  CuboidGeometry3D<T> cuboidGeometry(0,0,0,1, nx, ny, nz, noOfCuboids);

  cuboidGeometry.setPeriodicity( false, false, false );

  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  SuperGeometry3D<T> superGeometry( cuboidGeometry,loadBalancer);


  //OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;




  /*
   *  READING "F" STYLE
   *
   *
   */

	  T doubledummy;
	  ifstream numbersFile;
	  clout << "reading file ..." << std::endl;
	  numbersFile.open("LBgeometry_org.vtk");
	  if (!numbersFile.is_open())
	  {
		  clout << "reading file failed - program executed" << std::endl;
		  return 0;
	  }


	  int zz, yy, xx;
	  int locX, locY, locZ;
	  T noCuboid;
      clout << "bishierin" << std::endl;
	  for(xx=0; xx<nx; xx++){ //z loop
	    for(yy=0; yy<ny; yy++){ //y loop : y change faster than z
	      for(zz=0; zz<nz; zz++){ //x loop : x change faster than y
	    	  numbersFile >> doubledummy;
	    	 
	    	  noCuboid = superGeometry.getCuboidGeometry().get_iC(xx,yy,zz);


	    	  superGeometry.getCuboidGeometry().get(noCuboid).checkPoint(xx,yy,zz, locX, locY, locZ);
	    	  superGeometry.getBlockGeometry(noCuboid).get(locX,locY,locZ) = doubledummy;

	      }
	    }
	  }
	  clout << "reading file ... ok" << std::endl;


  	  //superGeometry.rename(1,2);
	  superGeometry.rename( 2,1, 1,1,1 ); // solid = 2

	  // inflow and outflow
	     Vector<T,3> origin(0,0,0);
	     Vector<T,3> normal(0,0,1);
	     T radiusBoundary = 100;
	     T eps = converter.getPhysDeltaX();
	     IndicatorCylinder3D<T> inflow( origin, normal, radiusBoundary, eps );
	     superGeometry.rename( 2,7,inflow );
	     superGeometry.rename( 7,3,1,1,0 );
	     superGeometry.rename(7,2);


	     origin[2] = nz-1;//-converter.getLatticeL();
	     normal[2] = -1;

	     IndicatorCylinder3D<T> outflow( origin, normal, radiusBoundary, eps );
	     superGeometry.rename( 2,8,outflow );
	     superGeometry.rename( 8,4,1,1,0 );
	     superGeometry.rename(8,2);




	  superGeometry.clean();
	  superGeometry.checkForErrors();

	  superGeometry.print();
	  clout << "Number of timesteps: " << converter.getLatticeTime( maxPhysT ) << std::endl;
	  clout << "press enter" << std::endl;
	  std::cin.ignore(10,'\n');




	  clout << "Prepare Geometry ... OK" << std::endl;
	  // inlet velocity (charU) depends on the number of inlet voxels:




	  // === 3rd Step: Prepare Lattice ===
	  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

	  //SmagorinskyBGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getOmega(),
	  //    instances::getBulkMomenta<T, DESCRIPTOR>(), 0.1, converter.getLatticeL(), converter.physTime() );
	  BGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getLatticeRelaxationFrequency(),
	                                          instances::getBulkMomenta<T, DESCRIPTOR>());

	



	  // choose between local and non-local boundary condition
	  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
	  createInterpBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryCondition );
	  // createLocalBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);

	  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
	  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );



	  prepareLattice( sLattice, bulkDynamics,
	                    sBoundaryCondition, sOffBoundaryCondition,
	                    superGeometry, converter );



	  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
	  timer.start();
	  util::ValueTracer<T> converge( converter.getLatticeTime(interval), epsilon );

	  for ( int iT = 0; iT <= converter.getLatticeTime( maxPhysT ); iT++ ) {
	


		  // === 5th Step: Definition of Initial and Boundary Conditions ===
		  setBoundaryValues( sLattice,converter, iT, superGeometry );

	      // === 6th Step: Collide and Stream Execution ===
	      sLattice.collideAndStream();

	      // === 7th Step: Computation and Output of the Results ===
	      getResults( sLattice, iT, superGeometry, converter);


	      if (sLattice.getStatistics().getAverageRho() != sLattice.getStatistics().getAverageRho())
	      {
	  		std::cout << "\n Rho diverges ... simulation killed \n" << std::endl;
	  		break;
	      }

	      if(converge.hasConverged()) {
	      	std::cout << "converged" << std::endl;
	      	break;
	      	}

	          converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
	          if(iT == converter.getLatticeTime( maxPhysT ) -1) clout << "ended normally" << endl;

	    }

	  timer.stop();
	  timer.printSummary();

}


