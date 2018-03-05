#include "SubsampleMain.h"
#include <stdlib.h>

SomeFunctions::SomeFunctions(void)
{
}

SomeFunctions::~SomeFunctions(void)
{
}

// =================================

// random number between 0 and 1
float SomeFunctions::randZeroOne(void)
{
	rand();
	return (float)rand()/RAND_MAX; // normalizes the random number
}

// decides if somethings happens with a given probability
bool SomeFunctions::isOccured(float probability)
{
	float randomNumber = randZeroOne();
	return ( randomNumber < probability );
}

void SomeFunctions::insertParameters(long &nX, long &nY, double &M, float &pF, bool &pF_x, short int &vd_type, short int &smpl_type, bool &ellip_mask, float& fully_sampled, bool flag_insertParam, short int &subsampleType)
{
	if(flag_insertParam)
	{
		bool rightInput = true;

		do
		{
			rightInput = true;

			cout << "=============================================================================="
			<< "\n Define subsampling strategy:"
			<< "\n 1: Poisson-Disc"
			<< "\n 2: Gaussian" << endl;

			cout << "Subsampling type [" << subsampleType << "]: ";
            subsampleType = getLineShort(subsampleType);
            cout << endl << endl;

            if(subsampleType == 1)
            {
                cout << "=============================================================================="
                << "\n--- Variable Density Poisson Disc Sampling ---"
                << "\n\n    parameters:"
                <<   "\n    -----------"
                <<   "\n      nX - the height of the image in [px] (between 50px and 512px)"
                <<   "\n      nY - the width of the image in [px] (between 50px and 512px)"
                <<   "\n      M  - the the downsamplingfactor of the sampling mask (between 2 and 15)"
                <<   "\n      pF - apply partial fourier/ESPReSSo to the sampling pattern"
                <<	 "\n           (1 (standard), 0.5, 0.625, 0.75 and 0.875)"
                <<   "\n      pF_x - orientation of partial fourier (1=x, 0=y)"
                <<   "\n      ellip_mask - elliptical mask switched on/off"
                <<   "\n      fully_sampled - percentag of fully sampled region"
                << "\n\n    available variable density types:"
                <<   "\n    ---------------------------------"
                << "\n        none = 1, centralPoint = 2, centralBlock = 3,"
                << "\n        centralEllipse = 4, fractionLoading = 5"
                << "\n\n    avaiable samplings types:"
                <<   "\n    -------------------------"
                << "\n        2D-PoissonSampling = 0, 3D-PoissonSampling = 1, \n"
                << "\n        3D-PureRandomSampling = 2\n"
                << "\n==============================================================================\n" << endl << endl << endl;
                cout << "Please enter the necessary parameters: " << endl << endl
				 << "nX [" << nX << "]       : ";
                nX = getLineLong(nX);
                cout << endl << "nY [" << nY << "]       : ";
                nY = getLineLong(nY);
                cout << endl << "M [" << M << "]          : ";
                M = getLineDouble(M);
                cout << endl << "pF [" << pF << "]         : ";
                pF = getLineFloat(pF);
                cout << endl << "pF_x [" << pF_x << "]       : ";
                pF_x = getLineBool(pF_x);
                cout << endl << "ellip_mask [" << ellip_mask << "] : ";
                ellip_mask = getLineBool(ellip_mask);
                cout << endl << "fully_sampled [" << fully_sampled << "]         : ";
                fully_sampled = getLineFloat(fully_sampled);
                cout << endl << "vd_type [" << vd_type << "]    : ";
                vd_type = getLineShort(vd_type);
                cout << endl << "smpl_type [" << smpl_type << "]  : ";
                smpl_type = getLineShort(smpl_type);
                cout << endl << endl;

            }
            else if(subsampleType == 2)
            {
                cout << "=============================================================================="
                << "\n--- Gaussian Sampling ---"
                << "\n\n    parameters:"
                <<   "\n    -----------"
                <<   "\n      nX - the height of the image in [px] (between 50px and 512px)"
                <<   "\n      nY - the width of the image in [px] (between 50px and 512px)"
                <<   "\n      M  - the the downsamplingfactor of the sampling mask (between 2 and 15)"
                <<   "\n      pF - apply partial fourier/ESPReSSo to the sampling pattern"
                <<	 "\n           (1 (standard), 0.5, 0.625, 0.75 and 0.875)"
                <<   "\n      pF_x - orientation of partial fourier (1=x, 0=y)"
                <<   "\n      fully_sampled - percentag of fully sampled region"
                << "\n==============================================================================\n" << endl << endl << endl;
                cout << "Please enter the necessary parameters: " << endl << endl
				 << "nX [" << nX << "]       : ";
                nX = getLineLong(nX);
                cout << endl << "nY [" << nY << "]       : ";
                nY = getLineLong(nY);
                cout << endl << "M [" << M << "]          : ";
                M = getLineDouble(M);
                cout << endl << "pF [" << pF << "]         : ";
                pF = getLineFloat(pF);
                cout << endl << "pF_x [" << pF_x << "]       : ";
                pF_x = getLineBool(pF_x);
                cout << endl << "fully_sampled [" << fully_sampled << "]         : ";
                fully_sampled = getLineFloat(fully_sampled);
                cout << endl << endl;

            }

			if( !((pF_x == 1) || (pF_x == 0)) || !(pF == 1 || pF == 0.5 || pF == 0.625 || pF == 0.75 || pF == 0.875) || nX <= 0 || nX > 1024 || nY <= 0 || nY > 1024
                    || M < 2 || M > 15 || vd_type < 1 || vd_type > 5 || smpl_type < 0 || smpl_type > 2 )
            {
                cout << "Wrong input!" << endl;
                rightInput = false;
            }
		} while(rightInput == false);
	}
	else
	{
		cout << endl << "--- Variable Density Poisson Disk Sampling ---" << endl
			 << "============================================" << endl << endl;
	}

	return;
}

long SomeFunctions::getLineLong(long defaultVal)
{
    string result;
    getline(cin, result);
    if(result.empty())
        return defaultVal;
    else
        return atol(result.c_str());
}

double SomeFunctions::getLineDouble(double defaultVal)
{
    string result;
    getline(cin, result);
    if(result.empty())
        return defaultVal;
    else
        return atof(result.c_str());
}

float SomeFunctions::getLineFloat(float defaultVal)
{
    string result;
    getline(cin, result);
    if(result.empty())
        return defaultVal;
    else
        return atof(result.c_str());
}

bool SomeFunctions::getLineBool(bool defaultVal)
{
    string result;
    getline(cin, result);
    if(result.empty())
        return defaultVal;
    else
        return result != "0";
}

short int SomeFunctions::getLineShort(short int defaultVal)
{
    string result;
    getline(cin, result);
    if(result.empty())
        return defaultVal;
    else
        return (short int)(atoi(result.c_str()));
}

void SomeFunctions::showParams(long nX, long nY, double M, float pF, bool pF_x, short int vd_type, short int smpl_type, bool ellip_mask, float fully_sampled, float deviation, short int body_region, short int subsampleType)
{
    cout << "==============================================================================" << endl;
    if(subsampleType == 1)
    {
        cout << "--- Generating Poisson-disc mask with ---" << endl;
    } else {
        cout << "--- Generating Gaussian mask with ---" << endl;
    }
    //<< "\n--- Starting mask generation with ---"
    cout << "\nsize: " << nX << " x " << nY
    << "\nacceleration: " << M << " (with deviation: " << deviation*M << " - " << (2-deviation)*M << ")"
    << "\nfully sampled: " << fully_sampled << endl;
    if(smpl_type == 0) {
        cout << "sampling type: " << smpl_type << " (2D)" << endl;
    } else {
        cout << "sampling type: " << smpl_type << " (3D)" << endl;
    }
    cout << "vd_type: " << vd_type << " | elliptical mask: " << ellip_mask
    << "\nPF/ESPReSSo: " << pF << " along " << pF_x
    << "\nBody region: " << body_region
	<< "\n==============================================================================" << endl;
}


