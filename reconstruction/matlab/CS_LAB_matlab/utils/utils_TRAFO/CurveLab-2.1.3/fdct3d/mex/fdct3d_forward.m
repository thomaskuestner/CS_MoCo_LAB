function C = fdct3d_forward(X)
  [m,n,p] = size(X);
  nbscales = floor(log2(min([m,n,p])))-2;
  nbdstz_coarse = 8;
  allcurvelets = 0;
  C = fdct3d_forward_mex(m,n,p,4,nbdstz_coarse,allcurvelets,X);
  
