function X = Inv_FineCurveCoeff(C);
% Inv_FineCurveCoeff: Reconstruct the last frequency subband from
%                           the finest scale curvelet coefficients
%  Usage:
%    X = Inv_FineCurveCoeff(C);
%  Inputs:
%    C   n by n matrix of fine scale curvelet coefficients
%  Outputs:
%    X   n by n matrix  
% See Also
%   FineCurveCoeff, Inv_Curvelet02Xform
%
% By Emmanuel Candes, 2003-2004
	  
	 X = fft2_mid0(C)/sqrt(prod(size(C)));
	 
	 
