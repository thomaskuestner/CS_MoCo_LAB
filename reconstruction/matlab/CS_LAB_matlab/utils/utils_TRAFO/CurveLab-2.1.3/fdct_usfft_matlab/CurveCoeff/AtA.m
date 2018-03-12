function Y = AtA(X,deep,IsImageReal);
% AtP: Gram Operator of SeparateAngle
%  Usage:
%     Y = AtA(X,deep)
%  Inputs:
%    X      n*n matrix (x,y)
%    deep   depth of angular localization
%  Outputs:
%    Y      n*n matrix (x,y)
%  Description:
%    Performs  A'A where A is the angular windowing 
% See Also
%   SeparateAngles, AtA_Toeplitz
%
% By Emmanuel Candes, 2003-2004

if nargin < 3,
  IsImageReal = 0;
end

        R  = SeparateAngles(X,deep,IsImageReal);
  	Y  = Adj_SeparateAngles(R,deep,IsImageReal);
      	
	
	
	
	 
	
	
 
