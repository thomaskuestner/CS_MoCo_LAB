function X = Inv_SeparateAngles(R,deep,IsImageReal,MaxIts,ErrTol,X0)
% Inv_SeparateAngles: 'Invert' Separate Angles by the method of
%                           Least Squares
%  Usage:
%    X = Inv_SeparateAngles(R,deep,MaxIts,ErrTol);
%  Inputs:
%    R       Array of smoothly localized Fourier samples
%    deep    Depth of the angular splitting
%    MaxIts  Maximum number of CG iterations. Default 10.
%    ErrTol  Error tolerance. Default 1.e-9.
%  Outputs:
%    X   n by n matrix  
%  Description
%    Performs the inverse angular separation. This is an
%    approximate inverse which uses a conjugate gradient solver on
%    the associated Gram system. The system is mildly preconditioned. 
% See Also
%   SeparateAngles, Inv_AtA_CG, 
%
% By Emmanuel Candes, 2003-2004

if nargin < 5,
	ErrTol = 1.e-9;
end
if nargin < 4,
	MaxIts = 10;
end
if nargin < 3,
  IsImageReal = 0;
end

n = size(R,3)*4;

W = SetScaleToZero(n);

P   = Adj_SeparateAngles(R,deep,IsImageReal).*W;
X   = Inv_AtA_CG(P,deep,MaxIts,ErrTol,X0);

