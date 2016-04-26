% dwt2D_plots
% DISPLAY 2D WAVELETS OF DWT2D.M

[af, sf] = farras;
J = 5;                      
L = 3*2^(J+1);
N = L/2^J;
x = zeros(L,3*L);
w = dwt2D(x,J,af);
w{J}{1}(N/2,N/2+0*N) = 1;
w{J}{2}(N/2,N/2+1*N) = 1;
w{J}{3}(N/2,N/2+2*N) = 1;
y = idwt2D(w,J,sf);

figure(1)
clf
imagesc(y);
axis equal
axis off
colormap(gray(128))

