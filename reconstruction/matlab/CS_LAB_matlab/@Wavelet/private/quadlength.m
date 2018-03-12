function [n,J] = quadlength(x)
% quadlength -- Find length and dyadic length of square matrix
%  Usage
%    [n,J] = quadlength(x)
%  Inputs
%    x   2-d image; size(n,n), n = 2^J (hopefully)
%  Outputs
%    n   length(x)
%    J   least power of two greater than n
%
%  Side Effects
%    A warning message is issue if n is not a power of 2,
%    or if x is not a square matrix.
%
	s = size(x);
	n = s(1);
	if s(2) ~= s(1),
		  disp('Warning in quadlength: nr != nc')
	end
	k = 1 ; J = 0; while k < n , k=2*k; J = 1+J ; end ;
	if k ~= n ,
		  disp('Warning in quadlength: n != 2^J')
	end

%
% Copyright (c) 1993. David L. Donoho
%     
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:40 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
