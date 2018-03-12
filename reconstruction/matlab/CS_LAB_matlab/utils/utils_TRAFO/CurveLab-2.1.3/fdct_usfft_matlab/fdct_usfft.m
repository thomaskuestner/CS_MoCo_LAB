function C = fdct_usfft(X,is_real,nscales)
% Fast Discrete Curvelet Transform via Unequispaced FFT's - Version 1.0
%
% Inputs
%   x           N by N pizel array (N dyadic) 
%
% Optional Inputs
%   is_real    Type of the transform
%                   0: complex-valued curvelets
%                   1: real-valued curvelets
%              [default set to 0]
%
%   nscales    number of scales including the coarsest wavelet level
%              [default set to log2(N) - 3, and we recommend
%              to take nscales <= default value. 
%   
% Outputs
%   C           Cell array of curvelet coefficients. 
%               C{j}{l}(k1,k2) is the coefficient at
%                   - scale j: integer, from coarsest to finest scale,
%                   - orientation l: labels the orientation of the
%                     curvelet (wedge), starting from the top-left
%                     wedge and increasing in a clockwise fashion
%                   - position k1,k2: both integers, size varies with j
%                     and l.
%
%               If is_real = 1, symmetric wedges around the origin
%               are paired as to return real-valued coefficients,
%               corresponding to  'cosine' or 'sine' curvelets
%
% See also ifdct_usfft.m, afdct_usfft.m
%
% By Emmanuel Candes, 2003--2004


  [m,n] = size(X);
  J = log2(n);
  
  if nargin < 3
    nscales = log2(n) - 3; 
  end
  
  if nargin < 2,
    is_real = 0;
  end
    
  L = J - nscales + 1; % Coarse scale
  scale = (L-1):(J-1);
  deep = floor(scale/2) + 1; 
  
  S = SeparateScales(X,L);
  
  C = cell(1,nscales);
  C{1} = cell(1,1);
  C{1}{1} =  CoarseCurveCoeff(S{1});
  for s = 2:nscales - 1
    C{s} = DetailCurveCoeff(S{s},deep(s),is_real); 
  end
  C{nscales} = cell(1,1);
  C{nscales}{1} = FineCurveCoeff(S{nscales});
  
  if is_real
    C = fdct_usfft_c2r(C);
  end
