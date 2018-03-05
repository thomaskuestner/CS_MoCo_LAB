#include "SubsampleMain.h"




VariableDensity::VariableDensity(void)
{
}

VariableDensity::VariableDensity(long nX, long nY, short int vd_type, float fully_sampled, bool ellip_mask, float p, float n, float iso_fac, double M)
{
	height = nX;
	width = nY;
	this->ellip_mask = ellip_mask;
	this->power = p;
	this->root = n;
	this->vd_type = vd_type;
	this->fully_sampled = fully_sampled;
	this->center[0] = long(floor(float(height)/2));
	this->center[1] = long(floor(float(width)/2));
	this->dRadius[0] = 0;
	this->dRadius[1] = 0;
	this->iso_fac = iso_fac;
	this->M = M;

	fraction = new float*[height];
	for (long i = 0; i < height; i++)
		fraction[i] = new float[width];

	for(long i=0; i<height; i++)
	{
		for(long j=0; j<width; j++)
		{
			fraction[i][j] = 1.0;
		}
	}
}

VariableDensity::~VariableDensity(void)
{
	for (int i = 0; i < height; i++)
		delete[] fraction[i];
	delete[] fraction;

	//delete [] dRadius;
	//delete [] center;
}

// =====================================================

float (**VariableDensity::getFraction(void))
{
	return fraction;
}
void VariableDensity::setFraction(long i, long j, float value)
{
	this->fraction[i][j] = value;
}
float VariableDensity::getFillvalue()
{
	return fillvalue;
}
void VariableDensity::setFillvalue(float fillvalue)
{
	this->fillvalue = fillvalue;
}
void VariableDensity::fsetRadius(double dRadius[])
{
	this->dRadius[0] = dRadius[0];
	this->dRadius[1] = dRadius[1];
}

// =====================================================

// choose the desired density option
long VariableDensity::genDensity()
{
	long nPointsMask;
	if (ellip_mask)
		nPointsMask = height * width - genEllipticalDen();
	else
		nPointsMask = height * width;

	//cout << "nPointsMask = "<< nPointsMask << endl;

	switch ( vd_type )
	{
	case 0: genLineSampling();
		break;
	case 1: genSingleDen();
		break;
	case 2: genCentralPointDen(nPointsMask);
		break;
	case 3: genCentralBlockDen(nPointsMask);
		break;
	case 4: genCentralEllipseDen(nPointsMask);
		break;
    case 5: genAprioriDen(nPointsMask);
        break;
	default:
		//cout << "Error: The desired Density Type is not avaiable." << endl;
		system("PAUSE");
		exit(1);
	}
	return nPointsMask;
}

long VariableDensity::genEllipticalDen()
{
	// generate the ellipse
	double radius_double[2] = { double(height) / 2, double(width) / 2 };
	double center_double[2] = { (double(height) / 2), (double(width) / 2)};

	long n = 0;
	for (long i = 0; i<height; i++)
	{
		for (long j = 0; j<width; j++)
		{
			if (1.0 < pow((i + 0.5 - center_double[0]) / radius_double[0], 2) + pow((j + 0.5 - center_double[1]) / radius_double[1], 2))
			{
				fraction[i][j] = -2;
				n++;
			}
		}
	}
	//cout << "number of -2s in ellipseMask: " << n << endl;
	return n;
}

// generate several densities
void VariableDensity::genLineSampling()
{
	center[1] = 0;

	// size of the fully sampled region block
	long fs = long(ceil(this->fully_sampled * float(height)));

	// create a line-density-map
	genFraction(center);

	fillvalue = fraction[center[0] - fs/2][0];
	if(fillvalue == 0)
		fillvalue = float(0.0001);

	// fill fraction with a marker at the positions of the fully sampled region
	for(long i = -fs/2; i <= fs/2; i++)
	{
		fraction[center[0]+i][0] = -1;
	}

	this->fully_sampled = float(fs +1)/height;
}

void VariableDensity::genSingleDen()
{
	fillvalue = 2;

	// fraction already initialised with ones in the constructor

	this->fully_sampled = 0;
}

void VariableDensity::genCentralPointDen(long nPointsMask)
{
	// generate density gradient
	genFraction(center);

	// compute the fillvalue with the mean of all 4 edges of the fully sampled region
	fillvalue = float(0.0001);

	// fill fraction with a marker at the positions of the fully sampled region
	fraction[center[0]][center[1]] = -1;

	this->fully_sampled = 1 / float(nPointsMask);
	//drawmask(height, width);
}

void VariableDensity::genCentralBlockDen(long nPointsMask)
{
	// size of the fully sampled region block
	long fs[2];
	fs[0] = long(ceil(sqrt(double(nPointsMask* fully_sampled*height) / double(width))));
	fs[1] = long(ceil(sqrt(double(nPointsMask* fully_sampled*width) / double(height))));

	// generate density gradient
	genFraction(center);

	// compute the fillvalue with the mean of all 4 edges of the fully sampled region
	if (fs[0] < fs[1])
		fillvalue = fraction[center[0] - fs[0] / 2][center[1]];
	else
		fillvalue = fraction[center[0]][center[1] - fs[1] / 2];

	// fill fraction and prob with a marker at the positions of the fully sampled region block
	for(int i = -fs[0]/2 ; i <= fs[0]/2 ; i++)
	{
		for(int j = -fs[1]/2 ; j<= fs[1]/2 ; j++)
		{
			fraction[center[0]+i][center[1]+j] = -1;
		}
	}
	this->fully_sampled = float((fs[0] + 1)*(fs[1] + 1))/float(nPointsMask);
	//drawmask(height, width);
}

void VariableDensity::genAprioriDen(long nPointsMask)
{
    // read in fraction in form of: -2 at center, -1 at fs region,  and distance normed to 1 at remaining pixels
    int i, j;
    ifstream in("fraction.txt");
    if (!in)
        {
            cout << "Cannot open fraction.txt.\n";
        }
    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            in >> fraction[i][j];
        }
    }
    in.close();


	fillvalue = float(0.0001);

    this->fully_sampled = fully_sampled;

	//drawmask(height, width);
}


void VariableDensity::genCentralEllipseDen(long nPointsMask)
{
	this->fully_sampled = determineDRadius(nPointsMask);

	long iX, iY;

	// size of the fully sampled region block
	double fs[2];
	if(dRadius[0] == 0 || dRadius[1] == 0)
	{
		fs[0] = double(ceil(fully_sampled * float(height)));
		fs[1] = double(ceil(fully_sampled * float(width)));
	}
	else
	{
		fs[0] = dRadius[0];
		fs[1] = dRadius[1];
	}

	// short axis of ellipse and fillvalue
	int nShortAxis;
	if (fs[0] < fs[1])
	{
		nShortAxis = int(ceil(fs[0]));
		fillvalue = fraction[center[0] - long(fs[0])][center[1]];
	}
	else
	{
		nShortAxis = int(ceil(fs[1]));
		fillvalue = fraction[center[0]][center[1] - long(fs[1])];
	}

	// generate density gradient
	genFraction(center);

	// generate the ellipse
	double center_double[2] = { (double(height) / 2), (double(width) / 2) };
	int nEllipsePoints = ellipse_grid_count(10*nShortAxis, fs, center_double);
	double *interior = ellipse_grid(10*nShortAxis, fs, center_double, nEllipsePoints);

	// fill fraction and prob with a marker at the positions of the fully sampled region
	long n = 0; //Ellipse Points in fully sampled region
	for(int i = 0; i<2*nEllipsePoints; i++)
	{
		// round the interior regions and mark them
		iX = long(interior[i]);
		iY = long(interior[i+1]);
		if(iX >= 0 && iX < height && iY >= 0 && iY < width)
		{
			if (fraction[iX][iY] != -1)
			{
				fraction[iX][iY] = -1;
				n++;
			}
		}
		i++;
	}

	delete[] interior;
	interior = 0;

	//cout << "number of -1s in mask: " << n << endl;

	this->fully_sampled = float(n)/float(nPointsMask);
	//drawmask(height, width);
}

// =====================================================

float VariableDensity::determineDRadius(long nPointsMask)
{
	float new_fully_sampled = fully_sampled;
	double fdeviation = 0.97;
	long **testFraction;
	testFraction = new long *[height];
	for( long i = 0 ; i < height ; i++ )
		testFraction[i] = new long[width];
	int lFixedLines = 0;
	double dFs = double(fully_sampled) * double(nPointsMask);
	dRadius[0] = double(ceil(fully_sampled * double(height)));
	dRadius[1] = double(ceil(fully_sampled * double(width)));

	for(long i = 0; i<height; i++)
	{
		for(long j = 0; j<width; j++)
		{
			testFraction[i][j] = 0;
		}
	}
	double dRadiusEst[2];
	dRadiusEst[0] = dRadius[0];
	dRadiusEst[1] = dRadius[1];
	double center_double[2] = { (double(height) / 2), (double(width) / 2)};
	int nShortAxis;
	int nEllipsePoints;
	double* interior;
	int iX, iY;

	if(vd_type == 4)
	{
		while(lFixedLines < long(fdeviation*dFs + 0.5))
		{
			lFixedLines = 0;
			for(long i = 0; i<height; i++)
			{
				for(long j = 0; j<width; j++)
				{
					testFraction[i][j] = 0;
				}
			}

			// short axis of ellipse
			if(dRadiusEst[0]<dRadiusEst[1])
				nShortAxis = int(ceil(dRadiusEst[0]));
			else
				nShortAxis = int(ceil(dRadiusEst[1]));

			nEllipsePoints = ellipse_grid_count(10*nShortAxis, dRadiusEst, center_double);
			interior = ellipse_grid(10*nShortAxis, dRadiusEst, center_double, nEllipsePoints);

			for(int i = 0; i<2*nEllipsePoints; i++)
			{
				iX = int(interior[i]);
				iY = int(interior[i+1]);
				if(iX >= 0 && iX < height && iY >= 0 && iY < width) {
					testFraction[iX][iY] = -1;
				}
				i++;
			}

			delete[] interior;
			interior = 0;

			for(int i = 0; i<height; i++)
			{
				for(int j = 0; j<width; j++)
				{
					if(testFraction[i][j] < 0)
						lFixedLines++;
				}
			}
			dRadius[0] = dRadiusEst[0];
			dRadius[1] = dRadiusEst[1];
			dRadiusEst[0] += 1; //double(ceil(double(height)/steps));
			dRadiusEst[1] += 1; //double(ceil(double(width)/steps));
		}

		//  delete testFraction
		for (int i = 0; i < height; i++)
			delete [] testFraction[i];
		delete[] testFraction;
		testFraction = 0;

		new_fully_sampled = float(lFixedLines)/float(nPointsMask);

		return new_fully_sampled;
	}
}


void VariableDensity::genFraction(long center[])
{
	// calculate the euclidean distance from each element of fraction to the center and insert it in fraction
		for (long i = 0; i<height; i++)
		{
			for (long j = 0; j<width; j++)
			{
				if (fraction[i][j] == 1)
				{
					float helper1;
					float helper2;
					helper1 = fabs(float(i) - float(center[0]));
					helper2 = fabs(float(j) - float(center[1]));
					fraction[i][j] = pow(pow(helper1, power) + pow(helper2, power), 1.0 / root);
				}
			}
		}


	// normalize all elements of fraction with its maximum value
	float fmax = findMaxOfFrac();
	for(long i = 0; i<height; i++)
	{
		for(long j = 0; j<width; j++)
		{
			if (fraction[i][j] != -2)
			{
				fraction[i][j] = fraction[i][j] / fmax;
			}
		}

	}
}

float VariableDensity::findMaxOfFrac(void)
{
	float fmax = 0;

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j<width; j++)
		{
			if (fraction[i][j] > fmax)
				fmax = fraction[i][j];
		}
	}
	return fmax;
}


double* VariableDensity::ellipse_grid ( int n, double r[2], double c[2], int ng )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_GRID generates grid points inside an ellipse.
//
//  Discussion:
//
//    The ellipse is specified as
//
//      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
//
//    The user supplies a number N.  There will be N+1 grid points along
//    the shorter axis.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R[2], the half axis lengths.
//
//    Input, double C[2], the center of the ellipse.
//
//    Input, int NG, the number of grid points inside the ellipse.
//
//    Output, double ELLIPSE_GRID[2*NG], the grid points.
//
{
  double h;
  int i;
  int j;
  int ni;
  int nj;
  int p;
  double x;
  double *xy;
  double y;

  xy = new double[2*ng];

  if ( r[0] < r[1] )
  {
    h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
    ni = n;
    nj = int( i4_ceiling ( r[1] / r[0] ) * ( double ) ( n ) );
  }
  else
  {
    h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
    nj = n;
    ni = int( i4_ceiling ( r[0] / r[1] ) * ( double ) ( n ) );
  }

  p = 0;

  for ( j = 0; j <= nj; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + ( double ) ( j ) * h;

	if (1.0 < pow((x - c[0]) / r[0], 2)
		+ pow((y - c[1]) / r[1], 2))
	{
		break;
	}

    xy[0+p*2] = x;
    xy[1+p*2] = y;
    p = p + 1;

    if ( 0 < j )
    {
      xy[0+p*2] = x;
      xy[1+p*2] = 2.0 * c[1] - y;
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + ( double ) ( i ) * h;

      if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 )
               + pow ( ( y - c[1] ) / r[1], 2 ) )
      {
        break;
      }

      xy[0+p*2] = x;
      xy[1+p*2] = y;
      p = p + 1;
      xy[0+p*2] = 2.0 * c[0] - x;
      xy[1+p*2] = y;
      p = p + 1;

      if ( 0 < j )
      {
        xy[0+p*2] = x;
        xy[1+p*2] = 2.0 * c[1] - y;
        p = p + 1;
        xy[0+p*2] = 2.0 * c[0] - x;
        xy[1+p*2] = 2.0 * c[1] - y;
        p = p + 1;
      }
    }
  }
  return xy;
}

int VariableDensity::ellipse_grid_count ( int n, double r[2], double c[2] )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPSE_GRID_COUNT counts the grid points inside an ellipse.
//
//  Discussion:
//
//    The ellipse is specified as
//
//      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
//
//    The user supplies a number N.  There will be N+1 grid points along
//    the shorter axis.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of subintervals.
//
//    Input, double R[2], the half axis lengths.
//
//    Input, double C[2], the center of the ellipse.
//
//    Output, int ELLIPSE_GRID)_COUNT, the number of grid points inside
//    the ellipse.
//
{
  double h;
  int i;
  int j;
  int ni;
  int nj;
  int p;
  double x;
  double y;

  if ( r[0] < r[1] )
  {
    h = 2.0 * r[0] / ( double ) ( 2 * n + 1 );
    ni = n;
    nj = int( i4_ceiling ( r[1] / r[0] ) * ( double ) ( n ) );
  }
  else
  {
    h = 2.0 * r[1] / ( double ) ( 2 * n + 1 );
    nj = n;
    ni = int( i4_ceiling ( r[0] / r[1] ) * ( double ) ( n ) );
  }

  p = 0;

  for ( j = 0; j <= nj; j++ )
  {
    i = 0;
    x = c[0];
    y = c[1] + ( double ) ( j ) * h;

    p = p + 1;

    if ( 0 < j )
    {
      p = p + 1;
    }

    for ( ; ; )
    {
      i = i + 1;
      x = c[0] + ( double ) ( i ) * h;

      if ( 1.0 < pow ( ( x - c[0] ) / r[0], 2 )
               + pow ( ( y - c[1] ) / r[1], 2 ) )
      {
        break;
      }

      p = p + 1;
      p = p + 1;

      if ( 0 < j )
      {
        p = p + 1;
        p = p + 1;
      }
    }
  }

  return p;
}

int VariableDensity::i4_ceiling ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    I4_CEILING rounds an R8 up to the next I4.
//
//  Example:
//
//    X        I4_CEILING(X)
//
//   -1.1      -1
//   -1.0      -1
//   -0.9       0
//   -0.1       0
//    0.0       0
//    0.1       1
//    0.9       1
//    1.0       1
//    1.1       2
//    2.9       3
//    3.0       3
//    3.14159   4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose ceiling is desired.
//
//    Output, int I4_CEILING, the ceiling of X.
//
{
  int value;

  value = ( int ) x;

  if ( value < x )
  {
    value = value + 1;
  }

  return value;
}

void VariableDensity::drawmask(long h, long w)
{
	ofstream fraction1;
	fraction1.open("fractionout.txt");

	//samplingPattern << endl << endl << "This is the sampling pattern: " << endl << endl;
	for (long i = 0; i < h; i++)
	{
		for (long j = 0; j < w; j++)
		{
			if ((i == center[0]) && (j == center[1]))
				fraction1 << "x";
			else if (fraction[i][j] == -1)
				fraction1 << "x";
			else if (fraction[i][j] == -2)
				fraction1 << "o";
			else
				fraction1 << "-";
		}
		fraction1 << endl;
	}
	fraction1.close();
}
