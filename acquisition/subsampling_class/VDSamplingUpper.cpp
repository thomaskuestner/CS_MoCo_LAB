#include "SubsampleMain.h"

VDSamplingUpper::VDSamplingUpper(void)
{
	long intOrigin[2] = {0,0};
	float floatOrigin[2] = {0,0};
	noPoint = Point(intOrigin, floatOrigin, false);
	anchor = 0;
	nElements = 0;
}

VDSamplingUpper::~VDSamplingUpper(void)
{
	for (int i = 0; i < height; i++)
		delete[] samplingMask[i];
	delete[] samplingMask;
	//samplingMask = 0;

	for (int i = 0; i < height; i++)
		delete[] grid[i];
	delete[] grid;
	//grid = 0;

	for (int i = 0; i < height; i++)
		delete[] grid2[i];
	delete[] grid2;
	//grid2 = 0;

    //delete [] pF_border;
}

VDSamplingUpper::VDSamplingUpper(bool flag_first, short int vd_type, short int smpl_type, long height, long width, double M, float pF_value, bool pF_x, float min_dist, int nPointsToTest, float deviation, long nPointsMask)

{
	// initialise parameters
	this->flag_first = flag_first;
	this->vd_type = vd_type;
	this->smpl_type = smpl_type;
	this->height = height;
	this->width = width;
	this->M = M;
	this->nPointsToTest = nPointsToTest;
	this->deviation = deviation;
	this->nPointsToCreate = (long)ceil(double(nPointsMask) / double(M));
	this->min_dist = min_dist;

	samplingMask = new int*[height];
	for (long i = 0; i < height; i++)
		samplingMask[i] = new int[width];

	grid = new Point*[height];
	for (long i = 0; i < height; i++)
		grid[i] = new Point[width];

	grid2 = new long*[height];
	for (long i = 0; i < height; i++)
		grid2[i] = new long[width];

	// initialise partial fourier matrix
	if(pF_x == true)
	{
		pF_border[0] = long(ceil((1-pF_value) * float(height)));
		pF_border[1] = 0;
	}
	else
	{
		pF_border[0] = 0;
		pF_border[1] = long(ceil((1-pF_value) * float(width)));
	}

	// initialise and resize the grids
	long intOrigin[2] = {0,0};
	float floatOrigin[2] = {0,0};
	noPoint = Point(intOrigin, floatOrigin, false);

	for(long i=0; i<height; i++){
		for(long j=0; j<width; j++){
			grid[i][j] = noPoint;
			grid2[i][j] = 0;
		}
	}
	nSamplingMaskPoints = 0;
	nSamplingMaskPointsHelper = 0;
	anchor = 0;
	nElements = 0;
	min_dist_status = 0;
}

// =================================

bool VDSamplingUpper::getFlag_first(void)
{
	return flag_first;
}
int VDSamplingUpper::getVd_type(void)
{
	return vd_type;
}
int VDSamplingUpper::getSmpl_type(void)
{
	return smpl_type;
}
double VDSamplingUpper::getM(void)
{
	return M;
}
long VDSamplingUpper::getNPointsToCreate(void)
{
	return nPointsToCreate;
}
int VDSamplingUpper::getNPointsToTest(void)
{
	return nPointsToTest;
}
float VDSamplingUpper::getMin_dist(void)
{
	return min_dist;
}
float VDSamplingUpper::getDeviation(void)
{
	return deviation;
}
long VDSamplingUpper::getHeight(void)
{
	return height;
}
long VDSamplingUpper::getWidth(void)
{
	return width;
}
int (**VDSamplingUpper::getSamplingMask(void))
{
	return samplingMask;
}
void VDSamplingUpper::getSamplingMask(int ***smplMask, long lPhase)
{
	for (int i = 0; i < height; i++)
	{
		for(int j = 0; j < width; j++)
		{
			smplMask[i][j][lPhase] = samplingMask[i][j];
		}
	}
}
long VDSamplingUpper::getNSamplingMaskPoints(void)
{
	return nSamplingMaskPoints;
}
long* VDSamplingUpper::getPFBorder(void)
{
	return pF_border;
}

void VDSamplingUpper::setMin_dist(float min_dist)
{
	this->min_dist = min_dist;
}
void VDSamplingUpper::setSamplingMask(int **samplingMask, long height, long width)
{
	for(long i=0; i<height; i++)
	{
		for(long j=0; j<width; j++)
		{
			this->samplingMask[i][j] = samplingMask[i][j];
		}
	}
}

// ====================================================================================================================================

void VDSamplingUpper::genMaskWithFullySampledRegion(VariableDensity *vd)
{
	//float fillvalue = vd->getFillvalue();
	// parameters to compute a new point
	Point newPoint;
	long newGridPos[2];
	float rCell;
	float lCell;
	float uCell;
	float dCell;

	// reininitialise the grids, lists and the sampling mask and fraction
	anchor = 0;
	nElements = 0;

	for (long k = 0; k<height; k++)
	{
		for (long j = 0; j<width; j++)
		{
			grid[k][j] = noPoint;
			samplingMask[k][j] = 0;
			grid2[k][j] = 0;

			if (vd->getFraction()[k][j] == -3)
				vd->setFraction(k, j, -1);
		}
	}


	nSamplingMaskPoints = 0;
	nSamplingMaskPointsHelper = 0;


	// fully sampled region
	////////////////////////

	for (long i = 0; i < height; i++)
	{
		for (long j = 0; j < width; j++)
		{
			if (vd->getFraction()[i][j] == -1)
			{

				// initialise the point
				newGridPos[0] = i;
				newGridPos[1] = j;
				newPoint.setGridPos(newGridPos);
				newPoint.pointLongToFloat(0.5, height, width);

				// save the point
				savePointToGrids(newPoint,vd);
				savePointInSmplMask(newPoint);

                if ((i >= pF_border[0]) && (j >= pF_border[1]))
				{
					// to avoid searching for a point, where it's not possible to find one-->
					// just save the point in the outer region of the fully sampled region:
					// four cells around the position of the newpoint in fraction
					rCell = vd->getFraction()[newGridPos[0]][min(newGridPos[1] + 1, width - 1)]; // right
					lCell = vd->getFraction()[newGridPos[0]][max(newGridPos[1] - 1, long(0))]; // left
					dCell = vd->getFraction()[min(newGridPos[0] + 1, height - 1)][newGridPos[1]]; // down
					uCell = vd->getFraction()[max(newGridPos[0] - 1, long(0))][newGridPos[1]]; // up
					if ((rCell != -1) || (lCell != -1) || (uCell != -1) || (dCell != -1))
						savePointInActList(newPoint);
				}
			}
		}
	}

	nSamplingMaskPointsHelper = nSamplingMaskPoints;

	cout << "nSamplingMaskPointsHelper = " << nSamplingMaskPointsHelper << endl;
	cout << "nElements = " << nElements << endl;

	for (int i = 0; i < nElements; i++)
	{
		Point actualPoint = LinkedList::showElement(i, anchor, nElements);
		vd->setFraction(actualPoint.getGridPos()[0], actualPoint.getGridPos()[1], -3.0);
	}
}

// decides which sampling will be used
short int VDSamplingUpper::genSamplingMask(VariableDensity *vd)
{
	switch(smpl_type)
	{
		case 0: return this->genLineSampling(vd); break;
		case 1: return this-> genPoissonSampling(vd); break;
		case 2: return this-> genRandomSampling(vd); break;
		default: //cout << "Error: The desired Sampling Type is not available." << endl;
				 system("PAUSE");
				 exit(1);
	}
	//return min_dist_status;
}


// ---This function generates the 2D Poisson Disk Sampling Pattern---
//////////////////////////////////////////////////////////////////////
short int VDSamplingUpper::genLineSampling(VariableDensity *vd)
{
	//setFractionWithVD(vd);
	float fillvalue = vd->getFillvalue();
	width = 1; // the sampling mask will be just a line in 2D Sampling

	// parameters to compute a new point
	Point newPoint;
	long newGridPos[2];
	Point active_points[2];
	active_points[0] = noPoint;
	active_points[1] = noPoint;

	// reininitialise the sampling mask
	for(long i=0; i<height; i++)
		samplingMask[i][0] = 0;

	nSamplingMaskPoints = 0;

	// 0 = undefined, 1=too small, 2=too large, 3=best approx, 4=done
	min_dist_status = 0;


	// fully sampled region
	////////////////////////

	if( vd_type != 1 )
	{
		for(long i=0; i<height; i++)
		{
			for(long j=0; j<width; j++)
			{
				if(vd->getFraction()[i][j] == -1)
				{
					if(i >= pF_border[0])
					{
						// initialise the point
						newGridPos[0] = i;
						newGridPos[1] = 0;
						newPoint.setGridPos(newGridPos);
						newPoint.pointLongToFloat(0.5, height, width); // or 0 instead of 0.5

						// save the point
						savePointInSmplMask(newPoint);

						// to avoid searching for a point, where it's no possible to find one-->
						// just save the point in the outer region of the fully sampled region:
						// four cells around the position of the newpoint in fraction
						float dCell = vd->getFraction()[newGridPos[0]+1][0]; // down
						float uCell = vd->getFraction()[newGridPos[0]-1][0]; // up
						if (dCell != -1)
							active_points[0] = newPoint;
						if (uCell != -1)
							active_points[1] = newPoint;
					}
				}
				vd->setFraction(active_points[0].getGridPos()[0], active_points[0].getGridPos()[1], fillvalue);
				vd->setFraction(active_points[1].getGridPos()[0], active_points[1].getGridPos()[1], fillvalue);
			}
		}

		// if just one fully sampled point is inserted
		if(active_points[1].getIsPoint() == 0)
			active_points[1] = active_points[0];
	}


	// first point
	///////////////

	if ( (active_points[0].getIsPoint() == 0) && (active_points[1].getIsPoint() == 0) )
	{
		bool found = false;
		do{
			// set the first point
			newGridPos[0] = height - 1 - int( floor(SomeFunctions::randZeroOne() * min_dist) );
			newGridPos[1] = 0;
			newPoint.setGridPos(newGridPos);
			newPoint.pointLongToFloat(-1, height, width);
			if( newGridPos[0] >= pF_border[0] )
				found = true;
		}while(!found);
		// save the first point
		savePointInSmplMask(newPoint);
		active_points[1] = newPoint;
	}


	// remaining points
	////////////////////

	// find points on the lower side
	while(active_points[0].getIsPoint() != 0)
	{
		newGridPos[0] = long(ceil( active_points[0].getGridPos()[0] + ( min_dist * vd->getFraction()[active_points[0].getGridPos()[0]][0] * (SomeFunctions::randZeroOne() + 1) ) ));
		newGridPos[1] = 0;
		newPoint.setGridPos(newGridPos);
		newPoint.pointLongToFloat(-1, height, width);

		// save the new point if it is inside the sampling line
		active_points[0] = newPoint;
		if( ( active_points[0].getRealPos()[0] > (height - 1) ) || ( active_points[0].getGridPos()[0] > (height - 1) ) )
			active_points[0] = noPoint;
		else
			savePointInSmplMask(newPoint);
	}

	// find points on the upper side
	while(active_points[1].getIsPoint() != 0)
	{
		newGridPos[0] = long(floor( active_points[1].getGridPos()[0] - ( min_dist * vd->getFraction()[active_points[1].getGridPos()[0]][0] * (SomeFunctions::randZeroOne() + 1) ) ));
		newGridPos[1] = 0;
		newPoint.setGridPos(newGridPos);
		newPoint.pointLongToFloat(-1, height, width);

		// save the new point if it is inside the sampling line
		active_points[1] = newPoint;
		if( ( active_points[1].getRealPos()[0] < 0 ) || ( active_points[1].getGridPos()[0] < 0 ) || ( active_points[1].getGridPos()[0] < pF_border[0] ) )
			active_points[1] = noPoint;
		else
			savePointInSmplMask(newPoint);
	}


	// Postproc step
	//////////////////

	if( nSamplingMaskPoints > ((2-deviation) * nPointsToCreate) )
	{
		min_dist_status = 1;
	}
	else if( nSamplingMaskPoints < (deviation * nPointsToCreate))
	{
		min_dist_status = 2;
	}
	else
	{
		min_dist_status = 3;
	}

	return min_dist_status;
}

// ---This function generates the standard 3D Poisson Disk Sampling Pattern---
///////////////////////////////////////////////////////////////////////////////
short int VDSamplingUpper::genPoissonSampling(VariableDensity *vd)
{
	float fillvalue = vd->getFillvalue();

	// indicate if a suitable mask was found and show the reason
	min_dist_status = 0;

	// parameters to compute a new point
	Point newPoint;
	long newGridPos[2];
	float newRealPos[2];
	long *newGridPosPtr;
	Point actualPoint;
	long *actGridPosPtr;
	float *actRealPosPtr;

	float radius;
	float angle;

	bool tooClose = true;

	int nextIndex = 0; // next index of active_list

	bool found = false;

	// reininitialise the grids, lists and the sampling mask
	anchor = 0;
	nElements = 0;
	for(long i=0; i<height; i++)
	{
		for(long j=0; j<width; j++)
		{
			if ((vd->getFraction()[i][j] != -1) && (vd->getFraction()[i][j] != -3))
			{
				grid[i][j] = noPoint;
				samplingMask[i][j] = 0;
				grid2[i][j] = 0;
			}
			if (vd->getFraction()[i][j] == -3) // border point
			{
				Point actualPoint = grid[i][j];
				filledCircle(actualPoint, vd);
				savePointInActList(actualPoint);
			}
		}
	}
	nSamplingMaskPoints = nSamplingMaskPointsHelper;

	// first point
	///////////////

	if (nElements == 0)
	{
		bool found_first_point = false;

		// create fix first sampling point to the center (not really center if height or width even!)
		if(flag_first)
		{
			// find first point
			newGridPos[0] = long(ceil(float(height)/2)-1);
			newGridPos[1] = long(ceil(float(width)/2)-1);

			// set first point, if partial fourier allows it
			if((newGridPos[0] >= pF_border[0]) && (newGridPos[1] >= pF_border[1]))
			{
				newPoint.setGridPos(newGridPos);
				newRealPos[0] = float(newGridPos[0])+0.5;
				newRealPos[1] = float(newGridPos[1])+0.5;
				newPoint.setRealPos(newRealPos);

				// save the first point
				savePoint(newPoint, vd);

				found_first_point = true;
			}
		}

		// create random first sampling point
		while(found_first_point == false)
		{
			// find the first point
			newRealPos[0] = float(height) * SomeFunctions::randZeroOne(); // was wenn randZeroOne = 1?
			newRealPos[1] = float(width) * SomeFunctions::randZeroOne();
			newPoint.setRealPos(newRealPos);
			newPoint.pointFloatToLong();
			newGridPosPtr = newPoint.getGridPos();

			// set first point, if partial fourier allows it
			if ((newGridPosPtr[0] >= pF_border[0]) && (newGridPosPtr[1] >= pF_border[1]) && vd->getFraction()[newGridPosPtr[0]][newGridPosPtr[1]] != -2)
			{
				// save the first point
				savePoint(newPoint, vd);

				found_first_point = true;
			}
		}
	}

	//cout << "nElements = " << nElements << endl;
	cout << ".";

	// remaining points
	////////////////////

	while( float(nSamplingMaskPoints)/float(nPointsToCreate) < (2.0 - deviation) )
	{

		found = false;

		if(nElements == 0)
		{
			min_dist_status = 2; // min_dist too large
			break;
		}
		else
		{
			// take a random element of the active list
			nextIndex = long(floor(0.5 + SomeFunctions::randZeroOne() * float(nElements-1)));
			actualPoint = LinkedList::showElement(nextIndex, anchor, nElements);
		}

		// randomly generate a new point around the active one allowed region
		// lies between minimal distance (for this voxel) and twice minimal distance
		for(int i = 0; i<nPointsToTest; i++)
		{
			actGridPosPtr = actualPoint.getGridPos();
			actRealPosPtr = actualPoint.getRealPos();

			// random radius depending on fraction
			if (vd->getFraction()[actGridPosPtr[0]][actGridPosPtr[1]] == -3)
			{
				radius = fillvalue * min_dist * (1 + SomeFunctions::randZeroOne());
			}
			else
				radius = vd->getFraction()[actGridPosPtr[0]][actGridPosPtr[1]] * min_dist * ( 1 + SomeFunctions::randZeroOne() );
			if (radius < 1) // avoids setting no point near a point with low values of fraction
				radius = 1;

			// random angle between 0..2pi
			angle = 2 * float(PI) * SomeFunctions::randZeroOne();

			// set position of the new point
			newRealPos[0] = actRealPosPtr[0] + radius * sin(angle);
			newRealPos[1] = actRealPosPtr[1] + radius * cos(angle);
			newPoint.setRealPos(newRealPos);

			// check if the point is in the grid and not inside the partial fourier region
			if ((newRealPos[0] >= pF_border[0]) && (newRealPos[1] >= pF_border[1]) && (newRealPos[0] < height) && (newRealPos[1] < width) && vd->getFraction()[long(newRealPos[0])][long(newRealPos[1])] != -2)
			{
				newPoint.pointFloatToLong();
				newGridPosPtr = newPoint.getGridPos();

				// check if there is already a point or a point disk
				if ( (grid[newGridPosPtr[0]][newGridPosPtr[1]].getIsPoint() == true) || (grid2[newGridPosPtr[0]][newGridPosPtr[1]] == 1))
					continue;

				// check if the newPoint is too close to other points
				checkNeighbourhood( newPoint, tooClose, vd );

				if ( !tooClose )
				{
					found = true;
					savePoint(newPoint,vd);
				}
			}
		}

		// remove actualPoint from active_list
		if ( /*!found &&*/ (nElements != 0) )
		{
				LinkedList::deleteElement(nextIndex, anchor, nElements);
		}
	}



	// Postproc step
	//////////////////

	// check if the number of generated points is in the deviation-range
	//if( (nElements == 0) && !(min_dist_status == 3) && ( float(nSamplingMaskPoints)/float(nPointsToCreate) >= deviation ) && ( float(nSamplingMaskPoints)/float(nPointsToCreate) <= (2-deviation) ) )
	//{
	//	min_dist_status = 3;
	//	cout << "...best approx found!!" << endl << endl;
	//}

	return min_dist_status;
}

// ---This function generates the randomSampling 3D Sampling Pattern---
///////////////////////////////////////////////////////////////////////////////
short int VDSamplingUpper::genRandomSampling(VariableDensity *vd)
{
	// generate a probability matrix which represents the probability
	// for placing a point in every single cell of the sampling mask
	float **prob = new float*[height];
		for( long i = 0 ; i < width; i++ )
			prob[i] = new float[512];
	if ( vd_type == 1 )
	{
		for(long i = 0; i <height; i++)
		{
			for(long j = 0; j <width-1; j++)
			{
				prob[i][j] = vd->getFraction()[i][j];
			}
		}
	}
	else
	{
		for(long i = 0; i <height; i++)
		{
			for(long j = 0; j <width-1; j++)
			{
				if(vd->getFraction()[i][j] == -1)
					prob[i][j] = 1;
				else
					prob[i][j] = 1-vd->getFraction()[i][j];
			}
		}
	}

	// indicate if a suitable mask was found and show the reason
	min_dist_status = 0;

	// parameters to compute a new point
	float probability;
	float probFactor = 3;
	long intPosition[2];

	// reininitialise the grids and lists
	for(long i=0; i<height; i++){
		for(long j=0; j<width; j++){
			samplingMask[i][j] = 0;
		}
	}
	nSamplingMaskPoints = 0;

	// fully sampled region
	////////////////////////

	if( vd_type != 1 )
	{
		for(long i=pF_border[0]; i<height; i++)
		{
			for(long j=pF_border[1]; j<width; j++)
			{
				if(vd->getFraction()[i][j] == -1)
				{
					// initialise the point
					intPosition[0] = i;
					intPosition[1] = j;

					// save the point
					savePointInSmplMask(intPosition);
				}
			}
		}
	}

	while ( nSamplingMaskPoints < nPointsToCreate )
	{
		// randomize a possible position for a point
		intPosition[0] = long(SomeFunctions::randZeroOne() * (height-1));
		intPosition[1] = long(SomeFunctions::randZeroOne() * (width-1));

		// check if a point already exists at the randomized position or if partial fourier declines the point
		if( (intPosition[0] < pF_border[0]) || (intPosition[1] < pF_border[1]) || (samplingMask[intPosition[0]][intPosition[1]] == 1) )
			continue;

		// decide if the point is used (depends on the probability matrix)
		probability = pow( prob[intPosition[0]][intPosition[1]] , probFactor );
		if( SomeFunctions::isOccured( probability ) )
			savePointInSmplMask(intPosition);
	}

	min_dist_status = 4;

	return min_dist_status;

}

// "visualize" the sampling pattern for the user in a .txt-file
void VDSamplingUpper::drawmask(long h, long w)
{
	ofstream samplingPattern;
	samplingPattern.open ("samplingPattern.txt");

	//samplingPattern << endl << endl << "This is the sampling pattern: " << endl << endl;
	for(long i=0; i<h; i++)
	{
		for(long j=0; j<w; j++)
		{
			if( (i == long(ceil(float(h)/2))-1) && (j == long(ceil(float(w)/2))-1) )
				samplingPattern << "0";
			else if( samplingMask[i][j] == 1 )
				samplingPattern << "x";
			else
				samplingPattern << "-";
		}
		samplingPattern << endl;
	}
	samplingPattern.close();
}


// checks, if any point is in the minimal Distance radius around a given points
void VDSamplingUpper::checkNeighbourhood(Point newPoint, bool &tooClose, VariableDensity *vd)
{
	tooClose = false;

	// parameters
	float *actRealPos;
	float *newRealPos;
	float minimalDistance;
	int gridMinDistance;
	long minCell[2];
	long maxCell[2];
	float c;
	long *gridPos = newPoint.getGridPos();

	minimalDistance = min_dist * vd->getFraction()[gridPos[0]][gridPos[1]];
	gridMinDistance = long(ceil(minimalDistance));

	// minimal cell position can be (1,1)
	minCell[0] = max(pF_border[0], gridPos[0]-gridMinDistance);
	minCell[1] = max(pF_border[1], gridPos[1]-gridMinDistance);

	// maximal cell index can be (height, width)
	maxCell[0] = min(height-1, gridPos[0]+gridMinDistance);
	maxCell[1] = min(width-1, gridPos[1]+gridMinDistance);

	// possible voxels to check which lie in the range of newPoint
	for( long i = minCell[0]; i <= maxCell[0]; i++)
	{
		for(long j = minCell[1]; j <= maxCell[1]; j++)
		{
			// there is a point in the actual cell
			if( samplingMask[i][j] == 1 )
			{
				actRealPos = grid[i][j].getRealPos();
				newRealPos = newPoint.getRealPos();
				c = sqrt( pow( actRealPos[0] - newRealPos[0], 2 ) + pow( actRealPos[1] - newRealPos[1], 2 ) );

				// new point too close to another one
				if( c < minimalDistance )
				{
					tooClose = true;
					return;
				}
			}
		}
	}

	return;
}

// draws a pseudo circle around the desired point into grid2
void VDSamplingUpper::filledCircle(Point thisPoint, VariableDensity *vd)
{
	// center of the circle
	long *center = thisPoint.getGridPos();

	// radius of the circle
	float r = min_dist * vd->getFraction()[center[0]][center[1]];

	int a;
	if (r < 0)
		a = 0;
	else // floored half edge length of the smallest square inside the circle
		a = int(floor(0.7071*r - 1));


	// fill grid2 with the 2a+1 x 2a+1 blocks
	for (int i = -a; i <= a; i++)
	{
		for (int j = -a; j <= a; j++)
		{
			if (center[0] + a >= 0 && center[0] + a <= height - 1 && center[1] + a >= 0 && center[1] + a <= width - 1)
				grid2[center[0] + a][center[1] + a] = 1;
		}
	}

	return;
}

// several functions to save a new point
void VDSamplingUpper::savePoint(Point p, VariableDensity *vd)
{
	long *gridPos = p.getGridPos();
	grid[gridPos[0]][gridPos[1]] = p;
	if (samplingMask[gridPos[0]][gridPos[1]] != 1)
	{
		filledCircle(p, vd); // fill grid2 with point and circle around it
		samplingMask[gridPos[0]][gridPos[1]] = 1;
		nSamplingMaskPoints++;
		LinkedList::insertElement(p, anchor, nElements);
	}
}
void VDSamplingUpper::savePointToGrids(Point p, VariableDensity *vd)
{
	long *gridPos = p.getGridPos();
	grid[gridPos[0]][gridPos[1]] = p;
	filledCircle(p, vd); // fill grid2 with point and circle around it
}
void VDSamplingUpper::savePointInActList(Point p)
{
	LinkedList::insertElement(p, anchor, nElements);
}
void VDSamplingUpper::savePointInSmplMask(Point p)
{
	long *gridPos = p.getGridPos();
	if (samplingMask[gridPos[0]][gridPos[1]] != 1)
	{
		samplingMask[gridPos[0]][gridPos[1]] = 1;
		nSamplingMaskPoints++;
	}
}
void VDSamplingUpper::savePointInSmplMask(long intPos[])
{
	if (samplingMask[intPos[0]][intPos[1]] != 1)
	{
		samplingMask[intPos[0]][intPos[1]] = 1;
		nSamplingMaskPoints++;
	}
}
