function y = windowing(x,L)
% This function computes L band sequence of a Meyer windowed signal.  
% Inputs:  x - signal to be decomposed
%          L - number of bandpass filters
% Outputs: y - windowed signal of x where second index indicates which 
%              bandpass index 
%
% Comments provided by Glenn R. Easley and Demetrio Labate. 
% Written by Wang-Q Lim on January 23, 2006.

N = length(x);
y = zeros(N,L);
T = N/L;
g = zeros(1,2*T);
for j = 1:2*T
    n = -1*T/2+j-1;
    g(j) = meyer_wind(n/T);
end

for j = 1:L
    index = 1;
    for k = -1*T/2:1.5*T-1
        in_sig=floor( mod(k+(j-1)*T,N))+1;
        y(in_sig,j) = g(index)*x(in_sig);
        index = index + 1;
    end
end
        
    

