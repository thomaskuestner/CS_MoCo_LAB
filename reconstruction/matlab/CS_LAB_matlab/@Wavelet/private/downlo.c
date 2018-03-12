#include "wavelab.h"

void downlo(x,n,lpf,m,y)
DOUBLE x[],lpf[],y[];
int   n,m;
{
           int n2, mlo, mhi, i, h, j; 
		   DOUBLE s;
           /*lowpass version */
           n2 = n/2;
		   mlo = m /2; 
		   mhi = n2-mlo;
		   if(2*mhi + (m-1) >= n) --mhi;
		   if(mhi < 0) mhi = -1;
           for( i= 0; i<=mhi; i++){
		   		s = 0.;
				for( h=0; h < m; h++)
					s += lpf[h]*x[2*i+h];
				y[i] = s;
			}
			
			/* fix up edge values */
			for( i= mhi+1; i<n2; i++){
		   		s = 0.;
				j = 2*i;
				for( h=0; h < m; h++){
					if(j >= n) j -= n;
					s += lpf[h]*x[j];
					j++;
				}
				y[i] = s;
			}
}
