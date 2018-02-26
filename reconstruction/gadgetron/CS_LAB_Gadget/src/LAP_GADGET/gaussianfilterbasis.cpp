#include "gaussianfilterbasis.h"

Gadgetron::GaussianFilterBasis::GaussianFilterBasis()
{
}

Gadgetron::GaussianFilterBasis::GaussianFilterBasis(int size_)
{
	size = ceil(size_);

	k.set_size(2*size+1, 2*size+1, 2*size+1);
	l.set_size(2*size+1, 2*size+1, 2*size+1);
	p.set_size(2*size+1, 2*size+1, 2*size+1);

	ColType k_ = linspace<ColType>(-size, size, fabs(-size-size)+1);
	MatType k__ = repmat(k_, 1, 2*size+1);

	for(int i = 0; i < k.n_slices; i++) {
		k.slice(i) = k__;
	}

	l = cs_lab_lap_arma::permuteSimple(k, 213);
	p = cs_lab_lap_arma::permuteSimple(k, 321);

	//Basis
	basis.set_size(pow(2*size+1,3), 4);

	CubeType b1 = cs_lab_lap_functions::g(k, size) % cs_lab_lap_functions::g(l, size) % cs_lab_lap_functions::g(p, size);
	b1 = b1 / accu(b1);
	CubeType b2 = cs_lab_lap_functions::g(k, size) % cs_lab_lap_functions::g(l, size) % cs_lab_lap_functions::g(p, size) % k;
	b2 = b2 / accu(b2 % k);
	CubeType b3 = cs_lab_lap_functions::g(k, size) % l % cs_lab_lap_functions::g(l, size) % cs_lab_lap_functions::g(p, size);
	b3 = b3 / accu(l % b3);
	CubeType b4 = cs_lab_lap_functions::g(k, size) % p % cs_lab_lap_functions::g(l, size) % cs_lab_lap_functions::g(p, size);
	b4 = b4 / accu(p % b4);

	for(int i = 0; i < b1.n_elem; i++) {
		basis(i) = b1(i);
	}

	for(int i = 0; i < b2.n_elem; i++) {
		basis(i+b1.n_elem) = b2(i);
	}

	for(int i = 0; i < b3.n_elem; i++) {
		basis(i+b1.n_elem+b2.n_elem) = b3(i);
	}

	for(int i = 0; i < b4.n_elem; i++) {
		basis(i+b1.n_elem+b2.n_elem+b3.n_elem) = b4(i);
	}

	//Seperabel filter
	G.set_size(2*size+1);
	Gi.set_size(2*size+1);
	Gd.set_size(2*size+1);
	Gdi.set_size(2*size+1);

	G = cs_lab_lap_functions::g(k_, size);

	Gd = k_ % G;
	G = G/arma::norm(G, 1);

	Gd = Gd/arma::norm(k_%Gd, 1);
	Gi = flipud(G);
	Gdi = flipud(Gd);
}

