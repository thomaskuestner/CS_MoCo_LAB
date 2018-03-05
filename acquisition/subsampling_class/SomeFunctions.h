#ifndef SOMEFUNCTIONS_H
#define SOMEFUNCTIONS_H

class SomeFunctions
{
public:
	SomeFunctions(void);
	~SomeFunctions(void);

// =================================

	// random number between 0 and 1
	static float randZeroOne(void);

	// decides if somethings happens with a given probability
	static bool isOccured(float probability);

	static void insertParameters(long &nX, long &nY, double &M, float &pF, bool &pF_x, short int &vd_type, short int &smpl_type, bool &ellip_mask, float& fully_sampled, bool flag_insertParam, short int &subsampleType);
	static void showParams(long nX, long nY, double M, float pF, bool pF_x, short int vd_type, short int smpl_type, bool ellip_mask, float fully_sampled, float deviation, short int body_region, short int subsampleType);

	static long getLineLong(long defaultVal);
	static double getLineDouble(double defaultVal);
	static float getLineFloat(float defaultVal);
	static bool getLineBool(bool defaultVal);
	static short int getLineShort(short int defaultVal);
};


#endif
