
% dualtree2D_plots
% DISPLAY 2D WAVELETS OF dualtree2D.M

J = 4;
L = 3*2^(J+1);
N = L/2^J;
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
x = zeros(2*L,3*L);
w = dualtree2D(x, J, Faf, af);
w{J}{1}{1}(N/2,N/2+0*N) = 1;
w{J}{1}{2}(N/2,N/2+1*N) = 1;
w{J}{1}{3}(N/2,N/2+2*N) = 1;
w{J}{2}{1}(N/2+N,N/2+0*N) = 1;
w{J}{2}{2}(N/2+N,N/2+1*N) = 1;
w{J}{2}{3}(N/2+N,N/2+2*N) = 1;
y = idualtree2D(w, J, Fsf, sf);

figure(1)
clf
imagesc(y);
axis equal
axis off
colormap(gray(128))

