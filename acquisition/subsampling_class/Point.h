#ifndef POINT_H
#define POINT_H


class Point
{
public:
	Point(void);
	Point(long gP[], float rP[], bool isP);
	~Point(void);

// =================================

	long* getGridPos(void);
	float* getRealPos(void);
	bool getIsPoint(void);

	void setGridPos(long gridPos[]);
	void setRealPos(float realPos[]);
	void setIsPoint(bool isPoint);

// =================================

	// "at = (0..1)" added to int-postion; (at<0)? then random number (0..1) is added
	static float pointLongToFloat(long intPos, float at, long height, long width);
	void pointLongToFloat(float at, long height, long width);

	static long pointFloatToLong(float floatPos);
	void pointFloatToLong(void);

	void printPoint(void); //prints to console
	
private:
	long gridPos[2];
	float realPos[2];
	bool isPoint;
};

#endif