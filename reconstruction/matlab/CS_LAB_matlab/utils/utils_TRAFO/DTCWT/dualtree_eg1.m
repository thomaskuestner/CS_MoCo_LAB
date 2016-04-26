J = 5;  % J: number of stages

% get filters
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

x = zeros(1,256);  % zero signal

% Compute dual-tree complex DWT of zero signal
w = dualtree(x, J, Faf, af); 
% Set a single (real) coefficient to 1
w{5}{1}(4) = 1;
% Compute the inverse transform 
y1 = idualtree(w, J, Fsf, sf);

% Compute dual-tree complex DWT of zero signal
w = dualtree(x, J, Faf, af); 
% Set a single (imaginary) coefficient to 1
w{5}{2}(4) = 1;
% Compute the inverse transform 
y2 = idualtree(w, J, Fsf, sf);

% Display real and imaginary parts and magnitude
n = [1:256]/256;
plot(n,y1,n,y2,n,sqrt(y1.^2+y2.^2))
title('COMPLEX 1D WAVELET') 
xlabel('t');
ylabel('\psi(t)');
