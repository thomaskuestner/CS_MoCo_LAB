/****************************************************************************************\
*	---Variable Density Subsampling---										 *
*	============================================										 *
*		available subsampling types:                                                     *
*		----------------------------    												 *
*       Poisson-Disc = 1, Gaussian = 2                                                   *
*																						 *
*		available variable density types:												 *
*		---------------------------------												 *
*  			single = 1, centralPoint = 2, centralBlock = 3								 *
*			centralEllipse = 4, localAdaptive = 5   				                     *
*																						 *
*		available Poisson-Disc sampling options (2D vs. 3D are chosen automatically):    *
*		---------------------------------------											 *
*   		2D-PoissonSampling = 0,	3D-PoissonSampling = 1,								 *
*			other 3D-Samplings: RandomSampling = 2										 *
*																						 *
*       ESPReSSo (Partial Fourier like compression):                                     *
*       --------------------------------------------                                     *
*           values = {0, 4/8, 5/8, 6/8, 7/8, 1}, direction: yPhase = false, zPhase = true*
*                                                                                        *
*       available generation of fraction: ((x-centerX)^p + (y-centerY)^p)^(1/n)   	     *
*       ---------------------------------                                                *
*           p = 0.5, 1, 2, 3, 4    n = 1, 2, p  (only LUTs for p = n = 1 or p = n = 2)   *
*																						 *
\****************************************************************************************/

// (c) copyright 2012 - 2016 under BSD license
// Thomas Kuestner, University of Tuebingen and University of Stuttgart, Germany
// (thomas.kuestner@{med.uni-tuebingen.de,iss.uni-stuttgart.de})

#include "SubsampleMain.h"

int main(int argc,char *argv[])
{
    // parameters
    long nX = 256;
    long nY = 128;
    short int vd_type = 4;
    short int body_region = 0;
    short int smpl_type = 1;
    float pF_value = 1;
    bool pF_x = false;
    bool ellip_mask = false;
    int lPhases = 1;
    double accels = 3;
    float p = 2;
    float n = 2;
    float iso_fac = 1.0;
    float fully_sampled = 0.065;
    short int subsampleType = 1;
    bool lWait = false;
    bool lPerform = true;

    if(argc > 1)
    {
        for(int i=1; i < argc; i++)
        {
            switch(i)
            {
                case 1:
                    if(std::string(argv[i]) == "--h" | std::string(argv[i]) == "-h" | std::string(argv[i]) == "--help" | std::string(argv[i]) == "-help")
                    {
                        cout << "============================================================================================================"
                        << "\n--- Subsampling Mask Generation ---"
                        << "\n============================================================================================================"
                        << "\n input parameters in ascending order (just numbers) with [default value] if left blank:"
                        << "\n no input parameters allows them to be entered in the command window"
                        << "\n\n parameter [defVal]\t\t\t description\t\t\t\t\t\t\t\t\t\t\t\t\t possible values"
                        << "\n------------------------------------------------------------------------------------------------------------"
                        << "\n nX ["<< nX <<"]:\t\t\t\t\t width of subsampling mask\t\t\t\t\t\t\t\t\t\t [16:1024]"
                        << "\n nY ["<< nY <<"]:\t\t\t\t\t height of subsampling mask\t\t\t\t\t\t\t\t\t\t [1:1024]"
                        << "\n accel ["<< accels <<"]:\t\t\t\t\t desired acceleration factor\t\t\t\t\t\t\t\t\t [2:15]"
                        << "\n subsampling type ["<< subsampleType <<"]:\t\t Poisson-Disc {1}, Gaussian {2}\t\t\t\t\t\t\t\t\t {1,2}"
                        << "\n fully sampled ["<< fully_sampled <<"]:\t\t percentage of fully sampled region\t\t\t\t\t\t\t\t [0:1]"
                        << "\n elliptical scanning ["<< ellip_mask <<"]:\t sample in elliptical scan region\t\t\t\t\t\t\t\t {0,1}"
                        << "\n ESPReSSo factor ["<< pF_value <<"]:\t\t ESPReSSo/Partial Fourier compactification factor\t\t\t\t [0.5:1]"
                        << "\n ESPReSSo direction ["<< pF_x <<"]:\t ESPReSSo/Partial Fourier direction along width {0}, height {1}\t {0,1}"
                        << "\n variable-density ["<<vd_type <<"]:\t\t variable-density options\t\t\t\t\t\t\t\t\t\t {1,2,3,4,5}"
                        << "\n \t\t\t\t\t\t\t none {1}"
                        << "\n \t\t\t\t\t\t\t central point {2}"
                        << "\n \t\t\t\t\t\t\t central block {3}"
                        << "\n \t\t\t\t\t\t\t central ellipse {4}"
                        << "\n \t\t\t\t\t\t\t local adaptive variable-density (needs external fraction.txt file) {5}"
                        << "\n Poisson Sampling ["<< smpl_type <<"]:\t\t Poisson-Disc sampling options\t\t\t\t\t\t\t\t\t {0,1,2}"
                        << "\n \t\t\t\t\t\t\t 2D Poisson-Disc {0} (chosen automatically if nY=1)"
                        << "\n \t\t\t\t\t\t\t 3D Poisson-Disc {1} (strictly obeying neighbourhood criterion)"
                        << "\n \t\t\t\t\t\t\t 3D pseudo Poisson-Disc {2} (approximately obeying neighbourhood criterion)"
                        << "\n power ["<< p <<"]:\t\t\t\t\t power of variable-density scaling\t\t\t\t\t\t\t\t [1:5]"
                        << "\n root ["<< n <<"]:\t\t\t\t\t root of variable-density scaling\t\t\t\t\t\t\t\t [1:5]"
                        << endl << endl;
                        lPerform = false;
                    } else {
                        nX = atol(argv[i]);
                    }
                    break;
                case 2: nY = atol(argv[i]); break;
                case 3: accels = double(atoi(argv[i])); break;
                case 4: subsampleType = atoi(argv[i]); break;
                case 5: fully_sampled = float(atof(argv[i])); break;
                case 6: std::stringstream(argv[i]) >> std::boolalpha >> ellip_mask; break;
                case 7: pF_value = float(atof(argv[i])); break;
                case 8: std::stringstream(argv[i]) >> std::boolalpha >> pF_x; break;
                case 9: vd_type = (short int)(atoi(argv[i])); break;
                case 10: smpl_type = (short int)(atoi(argv[i])); break;
                case 11: p = float(atof(argv[i])); break;
                case 12: n = float(atof(argv[i])); break;
                //case 13: body_region = (short int)(atoi(argv[i])); break;
                //case 14: iso_fac = float(atof(argv[i])); break;
                default: break;
            }
        }
    } else {
        SomeFunctions::insertParameters(nX, nY, accels, pF_value, pF_x, vd_type, smpl_type, ellip_mask, fully_sampled, true, subsampleType);
        lWait = true;
    }

    if(lPerform)
    {
        if(subsampleType == 1) // Poisson-disc
        {
            int ***samplingMask = startPD(nX, nY, accels, fully_sampled, pF_value, pF_x, lPhases, vd_type, smpl_type, ellip_mask, p, n, body_region, iso_fac);
        } else if(subsampleType == 2) // Gaussian
        {
            int ***samplingMask = startGaussian(nX, nY, accels, fully_sampled, pF_value, pF_x, long(lPhases));
        }
    }

    if(lWait)
        system("PAUSE");
    return EXIT_SUCCESS;
}


int ***startPD(long lLines, long lPartitions, double dAccel, float fully_sampled, float pF_val, bool pFx, int lPhases, short int v_type, short int s_type, bool ellipticalMask, float p, float n, short int body_part, float iso_fac)
{
    //////////////////////////////////////////
    //// Declarations and Initialisations ////
    //////////////////////////////////////////

    // seed timer for pseudo random number generation
    time_t seconds;
    time(&seconds);
    srand((unsigned) seconds);

    // flags
    bool flag_first = true;			  // fix first point on/off
    bool flag_autoTest = true;		  // automatic test if min_dist is too small or not on/off
    bool flag_insertParam = true;	  // insert parameters manually on/off

    // important parameters
    //vd_options vd_op;
    long nX = lLines;						// image height [px] y-direction
    long nY = lPartitions;					// image width [px] z-direction
    double M = dAccel;			  			// downsampling factor
    short int vd_type = v_type;		      // type of variable density map
    short int smpl_type = s_type;	      // type of sampling
    float pF_value = pF_val;		      // applies partial fourier to the poisson disk sampling
    bool pF_x = pFx;				      // orientation of the partial fourier
    bool ellip_mask = ellipticalMask;
    short int body_region = body_part;
    int ***samplingMask = 0;

    // options for variable density patterns
    //float fully_sampled = fs;		  // region of the image which is already fully sampled (in percent)

    unsigned int nPointsToTest = 20;  // determine how many points are generated/tested around one point
    float deviation = float(0.95);    // set maximal allowed deviation from undersampling factor

    float min_dist = M;				  // minimal distance around point in which no other point is allowed to be
    short int min_dist_status = 0;    // 0 = undefined, 1 = too small, 2 = too large, 3 = best approx, 4 = done

    // insert important parameters manually
    //SomeFunctions::insertParameters(nX, nY, M, pF_value, pF_x, vd_type, smpl_type, flag_insertParam);

    // avoid infinite calculation of min_dist if the number of fully sampled points
    // is equal or higher than the number of points to create
    if( (vd_type == 3 || vd_type == 4) && fully_sampled > float(1/M))
    {
        cout << "Error: number of fully_sampled-region-points is higher then the number of points to create!" << endl;
        system("PAUSE");
        exit(1);
    }
    if (body_region != 0 && vd_type == 1)
    {
        cout << "Error: Variable Density Type can't be 1 in an a priori mask"<< endl;
        system("PAUSE");
        exit(1);
    }

    long nY_helper = nY;				  // is '1' if 2D-Sampling is activated
    if(nY_helper == 1) // force 2D-Sampling
    {
        smpl_type = 0;
    }

    // 2D-Sampling - options
    if( smpl_type == 0 )
    {
        nY_helper = 1;
        ellip_mask = false;

        // there is just one variable density type for a 2D sampling mask
        //if( vd_type != 1 )
            vd_type = 0;
    }
    // 3D-Sampling - options
    else
    {
        if(vd_type == 0)
        {
            cout << "Variable Density Type can't be 0 in a 3D - Sampling Mask" << endl;
            system("PAUSE");
            exit(1);
        }
    }

    SomeFunctions::showParams(nX, nY, M, pF_value, pF_x, vd_type, smpl_type, ellip_mask, fully_sampled, deviation, body_region, 1);

    clock_t prgstart, prgende;		  // to measure the time

    // generate all density-informations for the sampling
    cout << "Generating Variable Density map..." << endl;
    VariableDensity* vd = new VariableDensity( nX, nY_helper, vd_type, fully_sampled, ellip_mask, p, n, iso_fac, M);

    if (vd_type < 6)
    {
        //if(smpl_type != 0) // 3D
            min_dist = Approx::findMinDistInLUT(nX, nY, M, fully_sampled, pF_value, vd_type, ellip_mask, p, iso_fac, smpl_type);
        //else
        //    min_dist = float(1/(2*M));
    }

    // necessary parameters for approximation of min_dist
    float step = 0;	 // change of min_dist
    //float range[2];  // range of values for min_dist
    /*float q;
    if (nX >= nY) q = float(nX) / float(nY);
    else q = float(nY) / float(nX);*/

    float* range;
    if(smpl_type != 0) {// 3D
        range = Approx::findRangeInLUT(nX, nY, M, fully_sampled, pF_value, vd_type, ellip_mask, p, iso_fac);
    } else {
        range = new float[2];
        range[0] = min_dist/100;
        range[1] = min_dist*6;
    }

    long nPointsMask = vd->genDensity();

    Approx *approx = new Approx(range, step);

    // generate Sampling object
    VDSamplingUpper *poiSamp = new VDSamplingUpper(flag_first, vd_type, smpl_type, nX, nY_helper, M, pF_value, pF_x, min_dist, nPointsToTest, deviation, nPointsMask);

	samplingMask = new int**[lLines];
	for (long i = 0; i < lLines; i++)
    {
        samplingMask[i] = new int*[lPartitions];
        for (long j = 0; j < lPartitions; j++)
            samplingMask[i][j] = new int[lPhases];
    }
    ///////////////////
    //// Main Loop ////
    ///////////////////

    cout << "\nProcessing Poisson Disc Sampling..." << endl;
    //cout << "... acceleration = " << M << endl << endl;
    //cout << "...starting with min_dist = " << min_dist << endl << endl;
    prgstart=clock();

    int lPhase;
    for (lPhase = 0; lPhase < lPhases; lPhase++)
    {
        int loopcounter = 0;
        bool failed = true;

        if (smpl_type == 1 && vd_type != 1)
            poiSamp->genMaskWithFullySampledRegion(vd);


        while (failed)
        {
            // generate the sampling mask
            min_dist_status = poiSamp->genSamplingMask(vd);

            // check the sampling mask and optimize min_dist if necessary
            failed = approx->checkMask(flag_autoTest, poiSamp, min_dist_status, vd);

            // too many iterations? --> stop loop
            loopcounter++;
            if (loopcounter == 150)
            {
                cout << "\nStopped main loop, number of iterations exceeded" << endl << endl;
                failed = false;
            }
        }


		// Postproc step
		//////////////////

		long nPoints = poiSamp->getNSamplingMaskPoints();					  // number of points in the sampling mask
		//int (**samplingMask)= poiSamp->getSamplingMask();
		poiSamp->getSamplingMask(samplingMask,lPhase);

		// draw the sampling mask
		poiSamp->drawmask(nX, nY);

		cout << "\n============================================" << endl << endl;
		cout << "Downsamplingfactor:   " << float(nPointsMask)/float(nPoints) << endl << endl;
		cout << "Sampled Points:       " << nPoints << endl << endl;
		cout << "Fully Sampled:        " << fully_sampled << endl << endl;
		prgende=clock();
		cout << "Time:                 "<< (float)(prgende-prgstart) / CLOCKS_PER_SEC <<" seconds"<< endl << endl;
		cout << "============================================" << endl << endl;
	}

    //delete poiSamp;
    poiSamp = 0;
    delete approx;
	delete vd;
    //delete &min_dist_status;

    return ( samplingMask );
}

int ***startGaussian(long lLines, long lPartitions, double dAccel, float fully_sampled, float pF_val, bool pFx, long lPhases)
{
    double		adSelectProb[MAX_PATTERN_LINES];
	double		dNormDist;
	double		dProbabilitySum = 0;
	double      dRadiusEst, dDistPart, dDistLine;
	long lLine , lPartition , lMeasurement, lFixedLines;
	long lEchoLine = long(lLines/2);
	long lEchoPartition = long(lPartitions/2);
	long m_lAcquiredLines	= long(double(lLines * lPartitions)/dAccel + 0.5);
    int ***samplingMask = 0;
    clock_t prgstart, prgende;

	// determine fixed lines
	if(lPartitions > 1) // 3D
    {
        // define local distance array and allocate memory
        double **dDistanceLocal = 0;
        dDistanceLocal = new double *[lLines];
        for(long i = 0; i < lLines; i++)
        {
            dDistanceLocal[i] = new double[lPartitions];
            for(long j = 0; j < lPartitions; j++)
            {
                dDistanceLocal[i][j] = 0;
            }
        }

		for(lPartition = 0; lPartition < lPartitions; lPartition++)
		{
			dDistPart = double(lPartition - lEchoPartition)/double(lPartitions);
			for(lLine = 0; lLine < lLines; lLine++)
			{
				//vDistance[lPartition*lLines + lLine] = 0.0; // initialize
				dDistLine = double(lLine - lEchoLine)/double(lLines);
				//vDistance[lPartition*lLines + lLine] = double(sqrt(pow(dDistPart, 2) + pow(dDistLine, 2))); // .push_back(double(sqrt(pow(dDistPart, 2) + pow(dDistLine, 2))));
				dDistanceLocal[lLine][lPartition] = double(sqrt(pow(dDistPart, 2) + pow(dDistLine, 2)));
            }
        }

		dRadiusEst = double(long(fully_sampled * 1000))/1000;
        while(lFixedLines < fully_sampled * lLines * lPartitions)
		{
			lFixedLines = 0;
			for(long i = 0; i < lLines; i++) //for(std::vector<int>::size_type i = 0; i != vDistance.size(); i++)
			{
				for(long j = 0; j < lPartitions; j++)
				{
					if(dDistanceLocal[i][j] <= dRadiusEst)
						lFixedLines++;
				}
			}
			dRadiusEst += 0.001;
        }

    } else {
		lFixedLines = long(fully_sampled * lLines);
		dRadiusEst = double(long(fully_sampled/2 * 1000))/1000;
    }


	samplingMask = new int**[lLines];
	for (long i = 0; i < lLines; i++)
    {
        samplingMask[i] = new int*[lPartitions];
        for (long j = 0; j < lPartitions; j++)
            samplingMask[i][j] = new int[lPhases];
    }

    short int smpl_type = 1;
    if(lPartitions == 1)
        smpl_type = 0;

    SomeFunctions::showParams(lLines, lPartitions, dAccel, pF_val, pFx, 0, smpl_type, false, fully_sampled, 1, 0, 2);
    // start mask generation process
    prgstart=clock();
	for(lLine = 0; lLine < lLines; lLine++) {
        dNormDist = fabs(double(lEchoLine - lLine) / double(lEchoLine + 10));
        adSelectProb[lLine] = pow(1 - dNormDist, 2);
        dProbabilitySum += adSelectProb[lLine];
    }

	long lInd = 0;
	long lCount;
	long alIndLine	[PDF_GENERATION_PREC];
    for(lLine = 0; lLine < lLines; lLine++) {
		lCount = long(adSelectProb[lLine]/dProbabilitySum*PDF_GENERATION_PREC + 0.5);
		for (long lLocalInd = 0; lLocalInd < lCount; lLocalInd++) {
			if (lLocalInd + lInd < PDF_GENERATION_PREC) alIndLine[lLocalInd + lInd] = lLine;
			}
        lInd += lCount;
    }
    long lTrueLUTLengthLine = (lInd<PDF_GENERATION_PREC)? lInd:PDF_GENERATION_PREC; // min()

    //. ---------------------------------------------------------------------------
	//. Create probability density function and LUT for partitions
	//. ---------------------------------------------------------------------------
	for(lPartition = 0; lPartition < lPartitions; lPartition++) {
		dNormDist = fabs(double(lEchoPartition - lPartition) / double(lEchoPartition + 5));
		adSelectProb[lPartition] = pow(1 - dNormDist, 2);
		dProbabilitySum += adSelectProb[lPartition];
	}

	lInd = 0;
	long alIndPart	[PDF_GENERATION_PREC];
	for(lPartition = 0; lPartition < lPartitions; lPartition++)
	{
		lCount = long(adSelectProb[lPartition]/dProbabilitySum*PDF_GENERATION_PREC + 0.5);
		for (long lLocalInd = 0; lLocalInd < lCount; lLocalInd++)
		{
			if (lLocalInd + lInd < PDF_GENERATION_PREC) alIndPart[lLocalInd + lInd] = lPartition;
		}
		lInd += lCount;
    }
	long lTrueLUTLengthPart = (lInd<PDF_GENERATION_PREC)? lInd:PDF_GENERATION_PREC; // min()

	//. ---------------------------------------------------------------------------
	//. Initialize the pattern by setting the fixed lines
	//. ---------------------------------------------------------------------------
	for (lPartition = 0; lPartition < lPartitions; lPartition++) {
		dDistPart = double(lPartition - lEchoPartition)/double(lPartitions);
		for(lLine = 0; lLine < lLines; lLine++) {
			dDistLine = double(lLine - lEchoLine)/double(lLines);
			if (sqrt(dDistPart*dDistPart + dDistLine*dDistLine) <= dRadiusEst) {
				for (lMeasurement = 0; lMeasurement < lPhases; lMeasurement++) samplingMask[lLine][lPartition][lMeasurement] = 1;
			} else {
				for (lMeasurement = 0; lMeasurement < lPhases; lMeasurement++) samplingMask[lLine][lPartition][lMeasurement] = 0;
			}
		}
	}

	//. ---------------------------------------------------------------------------
	//. Fill the sampling pattern line by line
	//. ---------------------------------------------------------------------------
	long lLineToActivate;
	long lPartToActivate;
	long lActiveLines;
	for (lMeasurement = 0; lMeasurement < lPhases; lMeasurement++) {
		lActiveLines = lFixedLines;
		while (lActiveLines < m_lAcquiredLines) {
			lInd = rand() % (lTrueLUTLengthPart);
			lPartToActivate = alIndPart[lInd];
			lInd = rand() % (lTrueLUTLengthLine);
			lLineToActivate = alIndLine[lInd];
			if (samplingMask[lLineToActivate][lPartToActivate][lMeasurement] == 0)
            {
				samplingMask[lLineToActivate][lPartToActivate][lMeasurement] = 1;
				lActiveLines++;
            }
        }
    }

    // print the sampling mask to file
    ofstream samplingPattern;
	samplingPattern.open ("samplingPattern.txt");

	//samplingPattern << endl << endl << "This is the sampling pattern: " << endl << endl;
    for(long k=0; k<lPhases; k++)
    {
        for(long i=0; i<lLines; i++)
        {
            for(long j=0; j<lPartitions; j++)
            {
                if( (i == long(ceil(float(lLines)/2))-1) && (j == long(ceil(float(lPartitions)/2))-1) )
                    samplingPattern << "0";
                else if( samplingMask[i][j][k] == 1 )
                    samplingPattern << "x";
                else
                    samplingPattern << "-";
            }
            samplingPattern << endl;
        }
        samplingPattern << endl << endl;
    }
	samplingPattern.close();

	cout << "\n============================================" << endl << endl;
    cout << "Downsamplingfactor:   " << float(lLines * lPartitions)/float(m_lAcquiredLines) << endl << endl;
    cout << "Sampled Points:       " << m_lAcquiredLines << endl << endl;
    cout << "Fully Sampled:        " << fully_sampled << endl << endl;
    prgende=clock();
    cout << "Time:                 "<< (float)(prgende-prgstart) / CLOCKS_PER_SEC <<" seconds"<< endl << endl;
    cout << "============================================" << endl << endl;

	return ( samplingMask );
}
