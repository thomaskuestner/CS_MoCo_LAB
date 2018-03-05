#ifndef VDSAMPLINGUPPER_H
#define VDSAMPLINGUPPER_H

class VDSamplingUpper
{
public:
	VDSamplingUpper(void);
	~VDSamplingUpper(void);
	VDSamplingUpper(bool flag_first, short int vd_type, short int smpl_type, long height, long width, double M, float pF_value, bool pF_x, float min_dist, int nPointsToTest, float deviation, long nPointsMask);

	// ==================================

	const double PI = 3.141592653589793;

	bool getFlag_first(void);
	int getVd_type(void);
	int getSmpl_type(void);
	double getM(void);
	long getNPointsToCreate(void);
	int getNPointsToTest(void);
	float getMin_dist(void);
	float getDeviation(void);
	long getHeight(void);
	long getWidth(void);
	int (**getSamplingMask(void));
	void getSamplingMask(int ***smplMask, long lPhase);
	long getNSamplingMaskPoints(void);
	long* getPFBorder(void);

	void setMin_dist(float m_d);
	void setSamplingMask(int **samplingMask, long height, long width);

	// =================================
    // prepare fully sampled region
	void genMaskWithFullySampledRegion(VariableDensity *vd);
	// generates sampling patterns (main method)
	short int genSamplingMask(VariableDensity *vd);
	// This function generates the 2D Poisson Disk Sampling Pattern
	short int genLineSampling(VariableDensity *vd);
	// This function generates the Poisson Disk Sampling Pattern
	short int genPoissonSampling(VariableDensity *vd);
	// This function generates the Random Sampling Pattern
	short int genRandomSampling(VariableDensity *vd);

	// ==================================

	// "visualize" the sampling pattern for the user in a .txt file
	void drawmask(long h, long w);

	// checks, if any point is in the minimal Distance radius around a given points
	void checkNeighbourhood(Point newPoint, bool &tooClose, VariableDensity *vd);

	// draws a pseudo circle around the desired point into grid2
	void filledCircle(Point thisPoint, VariableDensity *vd);

	// several functions to save a new point in grid and lists
	void savePoint(Point p, VariableDensity *vd);
	void savePointToGrids(Point p, VariableDensity *vd);
	void savePointInActList(Point p);
	void savePointInSmplMask(Point p);
	void savePointInSmplMask(long intPos[]);

	// ==================================

protected:

	long height;					  // image height [px] y-direction
	long width;					  // image width [px] z-direction
	double M;					  // downsampling factor

	short int vd_type;			  // type of variable density map
	short int smpl_type;		  // type of sampling

	bool flag_first;			  // fix first sampling point on/off

	long nPointsToCreate;		  // maximal points to create (global stopping condition)
	int nPointsToTest;			  // determine how many points are generated/tested around one point

	float min_dist;				  // minimal distance around point in which no other point is allowed to be
	float deviation;			  // set maximal allowed deviation from undersampling factor

	//// for all sampling algorithms
	int **samplingMask;	  // this is the sampling mask which will be computed
	long nSamplingMaskPoints;							   	  // counts the number of ones in the sampling mask
	long nSamplingMaskPointsHelper;							  // number of ones in sampling mask after filled with fully sampled region or after set first point
	short int min_dist_status;								  // shows if the minimal distance is too small, too large or ok
	long pF_border[2];										  // the points will be distributed untill this border (partial fourier)

	// for Poisson Sampling
	Point **grid; // to compute distance between points easily
	long **grid2;  // saves also circles around the points in which no other points can be
	//vector<Point> active_list;   // these points have to be processed
	LinkedList* anchor;								// anchor of the active list
	int nElements;									// number of elements in the active list
	Point noPoint;


};

#endif
