#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code;
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
//#include "../NMRShared.h"

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
   //lattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );
   //lattice.defineDynamics( superGeometry, 2 ,&bulkDynamics );
   lattice.defineDynamics( superGeometry,2,&instances::getNoDynamics<T,DESCRIPTOR>() );
   bc.addVelocityBoundary(superGeometry, 2, omega);





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

  lattice.defineRhoU( superGeometry,2,rhoF,uF );
    lattice.iniEquilibrium( superGeometry,2,rhoF,uF );


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
    sLattice.defineU( superGeometry, 3, uInlet);
    clout << "step=" << iT << "; maxVel=" << frac[0]*converter.getLatticeVelocity(converter.getCharPhysVelocity() ) << std::endl;


    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[2] = 2.25*frac[0]*converter.getLatticeVelocity(converter.getCharPhysVelocity() );

    T distance2Wall = converter.getPhysDeltaX()/2.;
    //RectanglePoiseuille3D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    //sLattice.defineU( superGeometry, 3, poiseuilleU );
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
  const int outputFlux = statIter;




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


//   if (iT%outputFlux == 0)
//   {
//
//	   std::vector<T> zPositions = {10., 75.};
//
//	   std::vector<int> materials = { 1, 3, 4 };
//	   for (auto& zPos : zPositions) {
//		   Vector<T, 3> center(nx/2,ny/2, zPos);
//		   Vector<T, 3> normal(0,0,1);
//		   T radius = 2./3.*nx;
//		   IndicatorCircle3D<T> slice( center, normal, radius );
//		   SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( sLattice, converter, superGeometry, slice, materials, BlockDataReductionMode::Discrete );
//		   int inputV[1] = {};
//		   T outputV[vFluxInflow.getTargetDim()] = {T()};
//		   vFluxInflow(outputV, inputV);
//		   cout << "Flux: " << outputV[0] << std::endl;
//		   vFluxInflow.print( "inflow","ml/s" );
//
//		   SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( sLattice, converter, superGeometry, slice, materials, BlockDataReductionMode::Discrete );
//		   int inputP[1] = {};
//		   T outputP[pFluxInflow.getTargetDim()] = {T()};
//		   pFluxInflow(outputP, inputP);
//		   cout << "FluxP: " << outputP[0] << std::endl;
//		   pFluxInflow.print( "inflow","N","mmHg" );
//
//	   }
//   }



}


int main( int argc, char* argv[] ) {



  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );


  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
      int {nx},                        // resolution: number of voxels per charPhysL
      (T)   0.607,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)   nx*deltaX,       // charPhysLength: reference length of simulation geometry
      (T)   0.0007,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)   1e-6, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)   998.0                       // physDensity: physical density in __kg / m^3__
    );

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

  CuboidGeometry3D<T> cuboidGeometry(0,0,0,converter.getPhysDeltaX(), nx, ny, nz, noOfCuboids);
//1
  cuboidGeometry.setPeriodicity( false, false, false );

  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  SuperGeometry3D<T> superGeometry( cuboidGeometry,loadBalancer);


  //OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;


//  // === Read MRI Data === //
//  // global analytical solution to optimize against
//  olb::AnalyticalF3D<T,T>* solution;
//
//  SuperIndicatorF3D<T>* objectiveDomain;
//    std::string dataName = "result.vti"; //opti_config["Application"]["MRIdata"].get<string>()+".vti";
//
//    BlockVTIreader3D<T,BaseType> readerData (dataName, "physVelocity");
//    BlockData3D<T,BaseType> data = readerData.getBlockData();
//    BlockDataF3D<T, BaseType> velocityF(data);
//    //BlockExtractIndicatorF3D<T> velocityLattice(velocityF, *objectiveDomain);
//
//    // Store the Analytical Functor of the velocity functor in the global variable solution
//    clout << "Interpolating over velocity functor..." << endl;
//    // Read Spacing_
//    XMLreader MRI_config(dataName);
//
//  //velocityLattice,readerData...
//    solution = new SpecialAnalyticalFfromBlockF3D<T,BaseType>(velocityF, readerData.getCuboid(), Vector<T,3> (2e-4, 2e-4, 2e-4));
//getLatticeU()
//    clout << "Interpolation finished." <<  endl;
//
//    solution->
//
//    std::cin.ignore(10, '\n');
  /*
   *  READING "F" STYLE
   *
   *
   */


	  T doubledummy;

	 	  clout << "reading file ..." << std::endl;
	 	  ifstream numbersFile;
	 	  if(singleton::mpi().getRank()==0){

	 	  numbersFile.open("LBgeometry.vtk");
	 	  if (!numbersFile.is_open())
	 	  {
	 		  clout << "reading file failed - program executed" << std::endl;
	 		  return 0;
	 	  }
	 	  }
	  int zz, yy, xx;
	  int locX, locY, locZ;
	  //T noCuboid;
	  T doubledummy2[1];
      clout << "bishierin" << std::endl;
	  for(xx=0; xx<nx; xx++){ //z loop
	    for(yy=0; yy<ny; yy++){ //y loop : y change faster than z
	      for(zz=0; zz<nz; zz++){ //x loop : x change faster than y
	    	  if(singleton::mpi().getRank()==0){
	    	    numbersFile >> doubledummy;

	    	    doubledummy2[0]=doubledummy;
	    	  }
	    	  clout << "doubledummy " << doubledummy2[0] << std::endl;
	    	  singleton::mpi().bCast(doubledummy2, 1);
	    	  singleton::mpi().barrier();
	    	  T physR[3];
	    	  physR[0]=xx*deltaX;
	    	  physR[1]=yy*deltaX;
	    	  physR[2]=zz*deltaX;
	    	  int latticeR[4];
	    	  if(superGeometry.getCuboidGeometry().getLatticeR(latticeR, physR )){
	    		  clout<<"error " <<latticeR[0]<<" "<<latticeR[1]<<" "<<latticeR[2]<<" "<<latticeR[3]<<std::endl;
	    		  clout<<"errorphys " <<physR[0]<<" "<<physR[1]<<" "<<physR[2]<<std::endl;
	    	      //superGeometry.set(latticeR[0],latticeR[1],latticeR[2],latticeR[3]) = doubledummy;
	    		if(loadBalancer.isLocal(latticeR[0])){
	    	      //superGeometry.getBlockGeometry(latticeR[0]).get(latticeR[1],latticeR[2],latticeR[3]) =1;
	    	      superGeometry.set(latticeR[0],latticeR[1],latticeR[2],latticeR[3]) = doubledummy2[0];
	    	      //superGeometry.getExtendedBlockGeometry(latticeR[0]).get(latticeR[1]+2,latticeR[2]+2,latticeR[3]+2) =1;
	    		}
	    	  }else {
	    		  clout<<"errorErroe" <<physR[0]<<" "<<physR[1]<<" "<<physR[2]<<std::endl;
	    	  }
	    	  //noCuboid = superGeometry.getCuboidGeometry().get_iC(xx,yy,zz);
	    	  //superGeometry.getCuboidGeometry().get(noCuboid).checkPoint(xx,yy,zz, locX, locY, locZ);
	    	  //superGeometry.getBlockGeometry(noCuboid).get(locX,locY,locZ) = doubledummy;
	      }
	    }
	  }
	  clout << "\n reading file ... ok" << std::endl;



//	  T doubledummy;
//	  ifstream numbersFile;
//	  clout << "reading file ..." << std::endl;
//	  numbersFile.open("LBgeometry.vtk");
//	  if (!numbersFile.is_open())
//	  {
//		  clout << "reading file failed - program executed" << std::endl;
//		  return 0;
//	  }
//
//
//	  int zz, yy, xx;
//	  int locX, locY, locZ;
//	  T noCuboid;
//      clout << "bishierin" << std::endl;
//	  for(xx=0; xx<nx; xx++){ //z loop
//	    for(yy=0; yy<ny; yy++){ //y loop : y change faster than z
//	      for(zz=0; zz<nz; zz++){ //x loop : x change faster than y
//	    	  numbersFile >> doubledummy;
//
//	    	  noCuboid = superGeometry.getCuboidGeometry().get_iC(xx,yy,zz);
//
//
//	    	  superGeometry.getCuboidGeometry().get(noCuboid).checkPoint(xx,yy,zz, locX, locY, locZ);
//	    	  superGeometry.getBlockGeometry(noCuboid).get(locX,locY,locZ) = doubledummy;
//
//	      }
//	    }
//	  }
//	  clout << "reading file ... ok" << std::endl;


  	  superGeometry.rename(1,2);
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


	     origin[2] = nz*deltaX-converter.getPhysDeltaX();
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


	  BGKdynamics<T, DESCRIPTOR> bulkDynamics(converter.getLatticeRelaxationFrequency(),
	                                          instances::getBulkMomenta<T, DESCRIPTOR>());

	



	  // choose between local and non-local boundary condition
	  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
	  createInterpBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryCondition );
	  //createLocalBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);

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


