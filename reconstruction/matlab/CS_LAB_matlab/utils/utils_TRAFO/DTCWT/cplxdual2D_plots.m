% cplxdual2D_plots
% DISPLAY 2D WAVELETS OF cplxdual2D.M

J = 4;
L = 3*2^(J+1);
N = L/2^J;
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
x = zeros(2*L,6*L);
w = cplxdual2D(x, J, Faf, af);
w{J}{1}{2}{2}(N/2,N/2+0*N) = 1;
w{J}{1}{1}{3}(N/2,N/2+1*N) = 1;
w{J}{1}{2}{1}(N/2,N/2+2*N) = 1;
w{J}{1}{1}{1}(N/2,N/2+3*N) = 1;
w{J}{1}{2}{3}(N/2,N/2+4*N) = 1;
w{J}{1}{1}{2}(N/2,N/2+5*N) = 1;
w{J}{2}{2}{2}(N/2+N,N/2+0*N) = 1;
w{J}{2}{1}{3}(N/2+N,N/2+1*N) = 1;
w{J}{2}{2}{1}(N/2+N,N/2+2*N) = 1;
w{J}{2}{1}{1}(N/2+N,N/2+3*N) = 1;
w{J}{2}{2}{3}(N/2+N,N/2+4*N) = 1;
w{J}{2}{1}{2}(N/2+N,N/2+5*N) = 1;
y = icplxdual2D(w, J, Fsf, sf);
y = [y; sqrt(y(1:L,:).^2+y(L+[1:L],:).^2)];
figure(1)
clf
imagesc(y);
title('2D Dual-Tree Complex Wavelets')
axis image
axis off
colormap(gray(128))
print -djpeg95 cplxdual2D_plots