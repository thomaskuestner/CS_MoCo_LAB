function X = Inv_CoarseCurveCoeff(C);
% Inv_CoarseCurveCoeff: Reconstruct the low-frequency subband from
%                           the coarsest scale curvelet coefficients
%  Usage:
%    X = Inv_CoarseCurveCoeff(C);
%  Inputs:
%    C   m by m matrix, m = 2*2^L
%  Outputs:
%    X   m by m matrix  
% See Also
%   CoaseCurveCoeff, Inv_Curvelet02Xform
%
% By Emmanuel Candes, 2003-2004

X = fft2_mid0(C)/sqrt(prod(size(C)));
	 
	 
	 
	 	 
