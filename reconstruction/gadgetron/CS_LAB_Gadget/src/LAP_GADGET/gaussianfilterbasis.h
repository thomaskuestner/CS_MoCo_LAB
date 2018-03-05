/*	
file name	: 	gaussianfilterbasis.h

author		: 	Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)

version		: 	1.0

date		: 	15.01.2018

description	: 	filter functions for LAP-based image registration
*/

#ifndef GAUSSIANFILTERBASIS_H
#define GAUSSIANFILTERBASIS_H
#include "mynamespaces.h"

using namespace arma;

namespace Gadgetron
{
	class GaussianFilterBasis
	{
	private:
		int size;
		CubeType k;
		CubeType l;
		CubeType p;

		MatType basis;

		//Seperable filters
		ColType G;
		ColType Gi;
		ColType Gd;
		ColType Gdi;

	public:
		GaussianFilterBasis();
		GaussianFilterBasis(int size);

		CubeType getK(){return k;}
		CubeType getL(){return l;}
		CubeType getP(){return p;}

		MatType getBasis(){return basis;}
		int getBasisWidth(){return basis.n_cols;}

		ColType getG(){return G;}
		ColType getGi(){return Gi;}
		ColType getGd(){return Gd;}
		ColType getGdi(){return Gdi;}
	};
}
#endif // GAUSSIANFILTERBASIS_H
