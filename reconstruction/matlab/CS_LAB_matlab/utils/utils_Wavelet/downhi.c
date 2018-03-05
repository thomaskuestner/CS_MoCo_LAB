void downhi(x,n,hpf,m,y)
double x[],hpf[],y[];
int   n,m;
{
           int n2, mlo, i, h, j; 
	   double s;

           /* highpass version */
           n2 = n/2;
           mlo = m/2-1; 
	   if(2*mlo+1 - (m-1) < 0) mlo++;
           for( i= mlo; i<n2; i++) {
                s = 0.;
                for( h=0; h < m; h++)
                     s += hpf[h]*x[2*i+1-h];
                y[i] = s;
           }
           if(mlo > n2) mlo = n2;
		/* fix up edge values */
           for( i= 0; i<mlo; i++) {
                s = 0.;
                j = 2*i+1;
                for( h=0; h < m; h++) {
                     if(j < 0) j += n;
                     s += hpf[h]*x[j];
                     --j;
                }
                y[i] = s;
           }
}
