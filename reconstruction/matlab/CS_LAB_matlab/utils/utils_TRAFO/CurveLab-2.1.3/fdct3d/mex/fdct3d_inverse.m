function X = fdct3d_inverse(C)
  [m,n,p] = size(C{end}{1});
  nbscales = floor(log2(min([m,n,p])))-2;
  nbdstz_coarse = 8;
  allcurvelets = 0;
  X = fdct3d_inverse_mex(m,n,p,4,nbdstz_coarse,allcurvelets,C);
 
 
