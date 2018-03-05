void uplo(x,n,lpf,m,y)
double x[],lpf[],y[];
int   n,m;
{
           int  meven, modd, i, h, j, mmax; 
		   double s;

           /*lowpass version */

			/* away from edges */
	       meven = (m+1)/2; modd = m/2;
           for( i= meven; i<n; i++){
		   		s = 0.;
				for( h=0; h < meven; h++)
					s += lpf[2*h]*x[i-h];
				y[2*i] = s;
				s = 0.;
				for( h=0; h < modd; h++)
					s += lpf[2*h+1]*x[i-h];
				y[2*i+1] = s;				
			}
			
			/* fix up edge values */
			mmax = meven;
			if(mmax > n) mmax = n;
			for( i= 0; i < mmax; i++){

				s = 0.;
				j = i;
				for( h=0; h < meven; h++){
					if(j < 0) j += n;
					s += lpf[2*h]*x[j];
					--j;
				}
				y[2*i] = s;

				s = 0.;
				j = i;
				for( h=0; h < modd; h++){
					if(j < 0) j += n;
					s += lpf[2*h+1]*x[j];
					--j;
				}
				y[2*i+1] = s;				
			}
}
			
          



			
          
