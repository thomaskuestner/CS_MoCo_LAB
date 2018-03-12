function z=conv2p(X,y)
% A function that performs periodic convolution 
%
% Input: X and y are the input arrays. X is assumed to be fft2 of x
%
% Output: z is the resultant convolution
%
% Written by Glenn R. Easley on Feb 2., 2006.
% Copyright 2006 by Glenn R. Easley. All Rights Reserved.
%

z=real(ifft2(X.*fft2(y)));
