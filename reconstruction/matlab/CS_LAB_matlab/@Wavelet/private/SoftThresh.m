function x = SoftThresh(y,t)
% SoftThresh -- Apply Soft Threshold 
%  Usage 
%    x = SoftThresh(y,t)
%  Inputs 
%    y     Noisy Data 
%    t     Threshold
%  Outputs 
%    x     sign(y)(|y|-t)_+
%
	res = (abs(y) - t);
	res = (res + abs(res))/2;
	x   = exp(i*angle(y)).*res;

%
% Copyright (c) 1993-5.  Jonathan Buckheit, David Donoho and Iain Johnstone
% modified by Michael Lustig for complex signals
%
    
    
 
 
%
%  Part of Wavelab Version 850
%  Built Tue Jan  3 13:20:39 EST 2006
%  This is Copyrighted Material
%  For Copying permissions see COPYING.m
%  Comments? e-mail wavelab@stat.stanford.edu 
