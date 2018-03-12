function X = afdct_usfft(C,isreal)

% afdct_usfft - Adjoint curvelet transform
%
% Input
%     C         Curvelet coefficients
%     isreal    Type of transform
%                   0: complex
%                   1: real
%
% Output
%     X         A double precision matrix
%

  [m,n] = size(C{end}{1});
  nbscales = floor(log2(min(m,n)))-3;
  nbangles_coarse = 16;
  allcurvelets = 0;
  
  if(isreal)
    C = fdct_usfft_c2r(C);
  end
  
  % call mex function
  X = afdct_usfft_mex(m,n,nbscales, nbangles_coarse, allcurvelets, C);
  

  
