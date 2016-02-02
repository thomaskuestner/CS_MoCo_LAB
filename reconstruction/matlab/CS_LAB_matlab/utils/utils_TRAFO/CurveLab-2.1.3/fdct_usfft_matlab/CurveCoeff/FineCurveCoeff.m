function C = FineCurveCoeff(X);
% FineCurveCoeff: Returns the tight frame coefficients at the finest scale
%  Usage:
%    C = FineCurveCoeff(X);
%  Inputs: 
%    X   n by n matrix. Typically, the input is the windowed Fourier
%        transform of an object. The window isolates the frequencies in
%        the last subband (finest scale).  
%  Outputs:
%    C   n by n matrix of wavelet coefficients  
%  Description:
%       Retrurns the tight frame coefficient table at the finest 
%       scale; these are "Mother wavelet coeffficients" associated
%       with frequency domain components near the dyadic corona
%       [-n, n]^2 \ [-n/4, n/4]^2. 
%
% By Emmanuel Candes, 2003-2004
  
  C = ifft2_mid0(X)*sqrt(prod(size(X)));
  
