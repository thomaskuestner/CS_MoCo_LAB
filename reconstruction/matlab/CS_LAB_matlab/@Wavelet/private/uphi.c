void uphi(x,n,hpf,m,y)
double x[],hpf[],y[];
int   n,m;
{
           int  meven, modd, i, h, j, mmin; 
		   double s;

		   /*hipass version */
	       meven = (m+1)/2; 
		   modd = m/2;

			/* away from edges */
           for( i= 0; i+meven<n; i++){
		   		s = 0.;
				for( h=0; h < meven; h++)
					s += hpf[2*h]*x[i+h];
				y[2*i+1] = s;
				s = 0.;
				for( h=0; h < modd; h++)
					s += hpf[2*h+1]*x[i+h];
				y[2*i] = s;				
			}
			
			/* fix up edge values */
			mmin = n-meven;
			if(mmin < 0) mmin = 0;
			for( i= mmin; i<n; i++){

				s = 0.;
				j = i;
				for( h=0; h < meven; h++){
					if(j >= n) j -= n;
					s += hpf[2*h]*x[j];
					j++;
				}
				y[2*i+1] = s;

				s = 0.;
				j = i;
				for( h=0; h < modd; h++){
					if(j >= n) j -= n;
					s += hpf[2*h+1]*x[j];
					j++;
				}
				y[2*i] = s;				
			}
}
