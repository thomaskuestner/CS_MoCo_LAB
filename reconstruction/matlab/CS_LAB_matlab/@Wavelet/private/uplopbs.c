#include <stdio.h>
void uplopbs(x,n,lpf,m,y)
double x[],lpf[],y[];
int   n,m;
{
  int  meven, modd, pos, i, h, j, k;
  double s;

  /*lopass version */
  meven = (m+1)/2; 
  modd = m/2;

  k = meven -1;
  j = (int) (k / 2);
  pos = k-2*j;
  k = j;

  for(i=0;i<n;i++){
    s = 0;
    j = k;
    for(h=0;h<(meven-pos);h++){
      while(j < 0) j += n;
      while(j >= n) j -= n;
      s += lpf[2*h+pos]*x[j];
      j--;
    }
    y[2*i] = s;

    s = 0;
    j = k+pos;
    for(h=0;h<(modd+pos);h++){
      while(j < 0) j += n;
      while(j >= n) j -= n;
      s += lpf[2*h+1-pos]*x[j];
      j--;
    }
    y[2*i+1] = s;

    k++;
  }
}
