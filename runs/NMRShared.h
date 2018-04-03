#ifndef NMRSHARED_H
#define NMRSHARED_H

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

void test_function(std::vector<double> vec)
{
	std::cout << "called the test function: " << vec[0] << std::endl;
	for (auto& zPos : vec) {
			   Vector<double, 3> center(nx/2,ny/2, zPos);
			   Vector<double, 3> normal(0,0,1);
			   T radius = 2./3.*nx;
			   IndicatorCircle3D<double> slice( center, normal, radius );
			   SuperPlaneIntegralFluxVelocity3D<double> vFluxInflow( sLattice, converter, superGeometry, slice, materials, BlockDataReductionMode::Discrete );
			   int inputV[1] = {};
			   T outputV[vFluxInflow.getTargetDim()] = {T()};
			   vFluxInflow(outputV, inputV);
			   cout << "Flux: " << outputV[0] << std::endl;
			   vFluxInflow.print( "inflow","ml/s" );

			   SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( sLattice, converter, superGeometry, slice, materials, BlockDataReductionMode::Discrete );
			   int inputP[1] = {};
			   T outputP[pFluxInflow.getTargetDim()] = {T()};
			   pFluxInflow(outputP, inputP);
			   cout << "FluxP: " << outputP[0] << std::endl;
			   pFluxInflow.print( "inflow","N","mmHg" );
}

#endif
