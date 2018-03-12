function X = ifdct_usfft(C,is_real)

% ifdct_usfft.m - Inverse Fast Discrete Curvelet Transform via
% Unequispaced FFTs - Version 1.0
%
% Inputs
%   C           Cell array containing the curvelet coefficient table (see
%               description in fdct_usfft.m)
%   is_real    Type of the transform
%                   0: complex-valued curvelets
%                   1: real-valued curvelets
%              [default set to 0]
% 
% Outputs
%   X           N-by-N matrix
% 
% Description
%    Performs the Least Squares inversion of the FDCT v USFFT. This
%    is an approximate inverse owing to the CG solver which is used
%    to unfold the angular partitioning
%
% See also fdct_usfft, afdct_usfft
%
% By Emmanuel Candes, 2003-2004

  [m,n] = size(C{end}{1});
  
  if is_real
    C = fdct_usfft_r2c(C);
  end
  
  nscales = length(C);
  S = cell(1,nscales);
  L = log2(n) - nscales + 1;
   
  S{1} = Inv_CoarseCurveCoeff(C{1}{1});
  
  for s = 2:nscales -1
    S{s} = Inv_DetailCurveCoeff(C{s},is_real);
      end
  S{nscales} = Inv_FineCurveCoeff(C{nscales}{1});
  
  X = Adj_SeparateScales(S,L);
  
  if (is_real)
  	 X = real(X);
  end
  
