function res = wavDenoise(im,lambda)

qmf = MakeONFilter('Daubechies',6);
M = FWT2_PO(im,1,qmf);
M_t = SoftThresh(M,lambda);
res = IWT2_PO(M_t,1,qmf);

