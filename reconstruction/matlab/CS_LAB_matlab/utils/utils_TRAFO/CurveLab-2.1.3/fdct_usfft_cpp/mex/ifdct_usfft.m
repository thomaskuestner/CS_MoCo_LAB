function X = ifdct_usfft(C, isreal)

% ifdct_usfft - Inverse curvelet transform
%
% Input
%     C         Curvelet coefficients
%     isreal    Type of transform
%                   0: complex
%                   1: real
%
% Output
%     X         Image
%

  [m,n] = size(C{end}{1});
  nbscales = floor(log2(min(m,n)))-3;
  nbangles_coarse = 16;
  allcurvelets = 0;
  
  if(isreal)
    C = fdct_usfft_r2c(C);
  end
  
  % call mex function
  X = ifdct_usfft_mex(m,n,nbscales, nbangles_coarse, allcurvelets, C);
  

  
