function [fxs,fys,fzs,nxs,nys,nzs] = fdct3d_param(C)
  [m,n,p] = size(C{end}{1});
  nbscales = floor(log2(min([m,n,p])))-2;
  nbdstz_coarse = 8;
  allcurvelets = 0;
  [fxs,fys,fzs,nxs,nys,nzs] = fdct3d_param_mex(m,n,p,nbscales,nbdstz_coarse,allcurvelets);
