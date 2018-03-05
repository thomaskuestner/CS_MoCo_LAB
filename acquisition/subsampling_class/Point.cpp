#include "SubsampleMain.h"

Point::Point(void)
{
	setIsPoint(true);
	this->gridPos[0] = 0;
    this->gridPos[1] = 0;
    this->realPos[0] = 0;
    this->realPos[1] = 0;
}

Point::Point(long gP[], float rP[], bool isP)
{
	setGridPos(gP);
	setRealPos(rP);
	setIsPoint(isP);
}

Point::~Point(void)
{
    //if(gridPos != 0)
    //    delete [] gridPos;
    //if(realPos != 0)
    //   delete [] realPos;
}

// =================================

long* Point::getGridPos(void)
{
	return gridPos;
}

float* Point::getRealPos(void)
{
	return realPos;
}

bool Point::getIsPoint(void)
{
	return isPoint;
}

void Point::setGridPos(long gridPos[])
{
	for(int i=0; i<=1; i++)
		this->gridPos[i] = gridPos[i];
}


void Point::setRealPos(float realPos[])
{
	for(int i=0; i<=1; i++)
		this->realPos[i] = realPos[i];
}

void Point::setIsPoint(bool isPoint)
{
	this->isPoint = isPoint;
}

// =================================

// "at = (0..1)" added to int-postion; (at<0)? then random number (0..1) is added
float Point::pointLongToFloat(long intPos, float at, long height, long width)
{
	float floatPos;


	if (at < 0) // randomize the position of the point
	{
        floatPos = (float(intPos) + SomeFunctions::randZeroOne());
	}
    else
	{
		floatPos = (intPos + at);
	}
	// avoid that the index exceeds the matrix dimensions
    if(floatPos > (height-1) || floatPos > (width-1))
		floatPos = float(intPos);

	return floatPos;
}

void Point::pointLongToFloat(float at, long height, long width)
{
	if (at < 0 || at >= 1) // randomize the position of the point
	{
		realPos[0] = (float(gridPos[0]) + SomeFunctions::randZeroOne());
		realPos[1] = (float(gridPos[1]) + SomeFunctions::randZeroOne());
	}
    else
	{
		realPos[0] = (gridPos[0] + at);
        realPos[1] = (gridPos[1] + at);
	}
	// avoid that the indices exceed the matrix dimensions
	if(realPos[0] > (height-1))
		realPos[0] = float((height-1));
	if(realPos[1] > (width-1))
		realPos[1] = float((width-1));
	return;
}

long Point::pointFloatToLong(float floatPos)
{
    long intPos = long(floatPos);
	return intPos;
}
void Point::pointFloatToLong(void)
{
	gridPos[0] = long(realPos[0]);
	gridPos[1] = long(realPos[1]);

	return;
}

void Point::printPoint(void)
{
	cout << "Point: grid position: ( " << gridPos[0] << " , " << gridPos[1] << " )"
		 << " and real position: ( " << realPos[0] << " , "  << realPos[1] << " )\n\n";
	return;
}
