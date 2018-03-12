function [index, window] = CoarseMeyerWindow(L,deg);
% CoarseMeyerWindow: Evaluates a lowpass Meyer window
%  Usage:
%    [index, window] = CoarseMeyerWindow(L,deg);
%  Inputs: 
%    L      Scale (integer)
%    deg    Degree of the polynomial
%  Outputs:
%    index  Location of the window on the time axis
%    window Low-frequency Meyer window

pio2 = pi/2;
dyadic_points = [1 2^L];

eps    = floor(dyadic_points(1)/3);
epsp   = floor(dyadic_points(2)/3);
lmidind = [ dyadic_points(1) : dyadic_points(2)-epsp ];
rmidind = [ dyadic_points(2)-epsp+1 : dyadic_points(2) ];
rghtind = [ dyadic_points(2)+1 : dyadic_points(2)+epsp+1 ];
farrghtind= [ dyadic_points(2)+epsp+2 : 2*dyadic_points(2)]; 

lmid = ones(1,length(lmidind));
rmid = cos(pio2*WindowMeyer((3/2)*((rmidind-1)/dyadic_points(2))-1,deg));
rght = cos(pio2*WindowMeyer((3/2)*((rghtind-1)/dyadic_points(2))-1,deg));
farrght =  zeros(1,length(farrghtind));

index = [lmidind rmidind rghtind farrghtind] -1;
window = [lmid rmid rght farrght];







