function [index, window] = FineMeyerWindow(dyadic_points,deg);
% FineMeyerWindow: Evaluates a highpass Meyer window
%  Usage:
%    [index, window] = FineMeyerWindow(L,deg);
%  Inputs: 
%    dyadic-points   Pair of the form 2^j, 2^(j+1)
%    deg    Degree of the polynomial
%  Outputs:
%    index  location of the window on the time axis; contains [2^j 2^(j+1)]
%    window highpass Meyer window

pio2 = pi/2;

eps    = floor(dyadic_points(1)/3);
epsp   = dyadic_points(1) - eps - 1;
farlftind = [ 1 : dyadic_points(1)-eps];
lftind     = [ dyadic_points(1)-eps+1 : dyadic_points(1)];
lmidind    = [ dyadic_points(1)+1 : dyadic_points(1)+eps+1];
rmidind    = [ dyadic_points(2)-epsp+1 : dyadic_points(2) ];

farlft = zeros(1,length(farlftind));
lft  = sin(pio2*WindowMeyer(3*((lftind-1)/dyadic_points(2))-1,deg));
lmid = sin(pio2*WindowMeyer(3*((lmidind-1)/dyadic_points(2))-1,deg));
rmid = ones(1,length(rmidind));

index = [farlftind lftind lmidind rmidind]; 
window = [farlft lft lmid rmid];



