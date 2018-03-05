#ifndef VARIABLEDENSITY_H
#define VARIABLEDENSITY_H

class VariableDensity
{
public:
	VariableDensity(void);
	VariableDensity(long nX, long nY, short int vd_type, float fully_sampled, bool ellip_mask, float p, float n, float iso_fac, double M);
	~VariableDensity(void);

// ==================================

	float **getFraction(void);
	float getFillvalue();
	void fsetRadius(double dRadius[]);
	void setFraction(long, long, float value);
	void setFillvalue(float fillvalue);
	float findMaxOfFrac(void);
// ==================================

	// chooses the desired density option
	long genDensity();
	long genEllipticalDen();

	// generate several densities
	void genLineSampling();
	void genSingleDen();
	void genCentralPointDen(long nPointsMask);
	void genCentralBlockDen(long nPointsMask);
	void genCentralEllipseDen(long nPointsMask);
	void genAprioriDen(long nPointsMask);


// ==================================

	float determineDRadius(long nPointsMask);

	// some other functions
	void genFraction(long center[]);
	double* ellipse_grid ( int n, double r[2], double c[2], int ng );
	int ellipse_grid_count ( int n, double r[2], double c[2] );
	int i4_ceiling ( double x );

	void drawmask(long h, long w);

private:

	float **fraction;		// an element of this matrix shows the distance to its center
	float fillvalue;		// value with which the fully sampled region will be filled later
	long height;			// height of fraction
	long width;				// width of fraction
	bool ellip_mask;
	long center[2];
	double dRadius[2];
	short int vd_type;		// variable density type
	float power;		    // generation of fraction: ((x-centerX)^power + (y-centerY)^power)^(1/root)
	float root;				// generation of fraction: ((x-centerX)^power + (y-centerY)^power)^(1/root)
	float iso_fac;          // anisotropic voxel dimension
	double M;

	// for centralBlock, centralEllipse
	float fully_sampled;	// region which is fully sampled in percent
};

#endif
