/*
file name	: 	LAP3D.h
author		: 	Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)
version		: 	1.0
date		: 	15.01.2018
description	: 	LAP-based image registration
*/

#ifndef LAP3D_H
#define LAP3D_H

#include <armadillo>
#include "cs_lab_lap_namespaces.h"

using namespace arma;

namespace Gadgetron
{
	class LAP3D {
	private:
		//Private member variables
		int m, n, p, numel;
		CubeType I1, I_input, I2, I2_shifted;
		ColType hGaussian, FilterSizes;
		field<CubeType> u_holder;
		CubeType mask[3];

		//Private member functions
		field<CubeType> estimateOpticalFlow3D(CubeType &I1_, CubeType &I2_, int K_);
		field<CubeType> estimateOpticalFlow3DKSpace(CubeType &I1_k_, CubeType &I2_k_, int sampleStepSize, int K_);

		ColType average(const ColType &I_, int K_);
		field<CubeType> cleanOF3D(field<CubeType> &u_est_);
		bool cleaningProcedure();

	public:
		LAP3D();
		LAP3D(const CubeType &I1_, const CubeType &I2_, int levelMin_, int levelMax_);

		field<CubeType> exec();
		
		void setMovingImage(const CubeType &I2_);
	};
}

#endif // LAP3D_H
