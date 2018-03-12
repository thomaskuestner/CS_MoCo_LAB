function y = meyer_wind(x)
% This function computes the Meyer window function of signal x
%
% Input   x - signal to be windowed
% 
% Output  y - windowed signal 
%
% Comments provided by Glenn R. Easley and Demetrio Labate. 
% Written by Wang-Q Lim on January 23, 2006.

if -1/3+1/2<x & x<1/3+1/2
    y = 1;
elseif (1/3+1/2<=x & x<=2/3+1/2) | (-2/3+1/2<=x & x<=1/3+1/2)
        w = 3*abs(x-1/2)-1;
        z = w^4*(35-84*w+70*w^2-20*w^3);
        y = cos(pi/2*(z))^2;
else
    y = 0;
end


    
        
        