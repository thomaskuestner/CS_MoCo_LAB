function C = fdct_usfft(X,isreal)

% fdct_usfft - Forward curvelet transform
%
% Input
%     X         Image
%     isreal    Type of transform
%                   0: complex
%                   1: real
% Output
%     C         Curvelet coefficients
%

  [m,n] = size(X);
  nbscales = floor(log2(min(m,n)))-3;
  nbangles_coarse = 16;
  allcurvelets = 0;
  
  %call mex function
  C = fdct_usfft_mex(m,n,nbscales, nbangles_coarse, allcurvelets, double(X));
  
  if(isreal)
    C = fdct_usfft_c2r(C);
  end
