void mirrorsymmfilt(lpf,hpf,length)
double *lpf, *hpf;
int    length;
{
    int i,k,isign;

    k = (length - 1) / 2;

    isign = 1;
    for(i=0; i<k;i++)
      isign *= -1;

    for(i=0; i < length; i++){
      *hpf++ = isign * *lpf++;
      isign *= -1;
    }
  }
