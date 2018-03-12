
#define DOUBLE double
void cpanalysis(double *x,double *cpt,double *temp,double *bm,double *bp,double *xt,double *y,int D,int n);
void fold(double *xc,double *xl,double *xr,double *bp,double *bm,int m,int n,double *y);
void makebell(int m,double *bp,double *bm);
void copydouble(DOUBLE *x,DOUBLE *y,int n);
void rshift(DOUBLE *x,DOUBLE *y,int n);
void half(DOUBLE *x,int n);
void lshift(DOUBLE *x,DOUBLE *y,int n);
void copyadd(DOUBLE *x,DOUBLE *y,int n);
void adddouble(DOUBLE *x,DOUBLE *y,int n, DOUBLE *z);
void unpackdouble(DOUBLE *x,int n,int nc,int k,DOUBLE *y);
void packdouble(DOUBLE *x,int n,int nc,int k,DOUBLE *y);
void copyreverse(DOUBLE *x,DOUBLE *y,int n);
void changesign(double *x,int n);
void dctiv(double *x,double *y,double *t,int n);
void four1(double *data,unsigned int nn,int isign);
void realft(double *data,unsigned int n, int isign);
void wpd(double *sig,int nr,int Dee,double *hpf,double *lpf,int lenfil,double *wc,double *temp);
void downhi(double *x,int n,double *hpf,int m,double *y);
void downlo(double *x,int n,double *lpf,int m,double *y);
void mirrorfilt(double *lpf,double *hpf,int length);

void fillzeros(double *x,int n);
void edgeadd(int which,double *x,double *y,int n,int m);
void edgefold(int which,double *x,double *y,int n,int m,double *bp,double *bm);
void unfold(double *S,double *xc,double *xl,double *xr,int n,int m,double *bp,double *bm);
void edgeunfold(int which,double *x,int n,int m,double *bp,double *bm);

#define FLOAT double
void hartley(
  FLOAT  *x,
  int     m,
  FLOAT  *tab);

void dct(
  FLOAT *x,
  FLOAT *y,
  int    m,
  FLOAT *tab);
  
void downhipbs(double *x,int n,double *hpf,int m,double *y);
void downlopbs(double *x,int n,double *lpf,int m,double *y);
void dst(double *x,double *y,int m,double *tab);

void idct(
  FLOAT *x,
  FLOAT *y,
  int    m,
  FLOAT *tab);

void idst(double *x,double *y,int m,double *tab);
void matinv(double *a,unsigned int n);
void maiseg(double *eta0,double *SegFiltLeft,double *SegFiltRight);
void matmpy(double *a,int n,int m,double *b,int l,double *c);
void mirrorsymmfilt(double *lpf,double *hpf,int length);
void uphi(double *x,int n,double *hpf,int m,double *y);
void uplo(double *x,int n,double *lpf,int m,double *y);
void uphipbs(double *x,int n,double *hpf,int m,double *y);
void uplopbs(double *x,int n,double *lpf,int m,double *y);



