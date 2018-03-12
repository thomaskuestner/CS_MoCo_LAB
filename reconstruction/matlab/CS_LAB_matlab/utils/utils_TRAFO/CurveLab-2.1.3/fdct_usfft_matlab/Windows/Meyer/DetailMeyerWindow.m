function [index, window] = DetailMeyerWindow(dyadic_points,deg);
% DetailMeyerWindow: Evaluates a bandpass Meyer window
%  Usage:
%    [index, window] = DetailMeyerWindow(L,deg);
%  Inputs: 
%    dyadic-points   Pair of the form 2^j, 2^(j+1)
%    deg    Degree of the polynomial
%  Outputs:
%    index  location of the window on the time axis; contains [2^j 2^(j+1)]
%    window bandpass Meyer window


pio2 = pi/2;

eps    = floor(dyadic_points(1)/3);
epsp   = dyadic_points(1) - eps - 1;
lftind  = [ dyadic_points(1)-eps+1 : dyadic_points(1) ];
lmidind = [ dyadic_points(1)+1 : dyadic_points(1)+eps+1 ];
rmidind = [ dyadic_points(2)-epsp+1 : dyadic_points(2) ];
rghtind = [ dyadic_points(2)+1 : dyadic_points(2)+epsp+1 ];

lft  = sin(pio2*WindowMeyer(3*((lftind-1)/dyadic_points(2))-1,deg));
lmid = sin(pio2*WindowMeyer(3*((lmidind-1)/dyadic_points(2))-1,deg));
rmid = cos(pio2*WindowMeyer((3/2)*((rmidind-1)/dyadic_points(2))-1,deg));
rght = cos(pio2*WindowMeyer((3/2)*((rghtind-1)/dyadic_points(2))-1,deg));

index = [lftind lmidind rmidind rghtind] -1;
window = [lft lmid rmid rght];



