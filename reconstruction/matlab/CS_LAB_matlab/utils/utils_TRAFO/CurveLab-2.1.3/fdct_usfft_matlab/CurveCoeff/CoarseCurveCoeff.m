function C = CoarseCurveCoeff(X);
% CoarseCurveCoeff: Returns the tight frame coefficients 
% at the coarsest scale;
%  Usage:
%      C = CoarseCurveCoeff(X);
%  Inputs:
%    X      m by m matrix, m = 2*2^L. Typically, the input is
%           obtained after multiplying the Fourier transform of an object
%           with a low frequency Meyer window   
%  Outputs:
%    C      m by m matrix  
%  Description:
%       Returns the tight frame coefficient table at the
%       coarsest scale; these are "Father wavelet coeffficients," 
%       with frequency domain components near the dyadic
%       low-frequency square Q = [-2^(L-1), 2^(L-1)]^2.
% See Also
%   Curvelet02Xform, Inv_CoarseCurveCoeff, SeparateScales
%
% By Emmanuel Candes, 2003-2004

	 C = ifft2_mid0(X)*sqrt(prod(size(X)));
     
