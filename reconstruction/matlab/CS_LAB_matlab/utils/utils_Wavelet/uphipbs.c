void uphipbs(x,n,hpf,m,y)
double x[],hpf[],y[];
int   n,m;
{
  int  meven, modd, pos, i, h, j, k;
  double s;

  /*hipass version */
  meven = (m+1)/2; 
  modd = m/2;

  k = meven + 2*n - m;
  j = (int) (k / 2);
  pos = k-2*j;
  k = j;

  for(i=0;i<n;i++){
    s = 0;
    j = k + pos;
    for(h=0;h<meven-pos;h++){
      while(j >= n) j -= n;
      while(j < 0) j += n;
      s += hpf[2*h+pos]*x[j];
      j++;
    }
    y[2*i+1] = s;

    s = 0;
    j = k;
    for(h=0;h<modd+pos;h++){
      while(j >= n) j -= n;
      while(j < 0) j += n;
      s += hpf[2*h+1-pos]*x[j];
      j++;
    }
    y[2*i] = s;
    
    k++;
  }
}
