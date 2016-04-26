% dwt3D_plots
% DISPLAY 3D WAVELETS OF DWT3D.M

[af, sf] = farras;
J = 4;
L = 3*2^(J+1);
N = L/2^J;
x = zeros(L,L,L);
w = dwt3D(x,J,af);
w{J}{7}(N/2,N/2,N/2) = 1;
y = idwt3D(w,J,sf);
figure(1)
clf
v = 1:L;
S = 0.017;
p1 = patch(isosurface(v,v,v,y,S));
isonormals(v,v,v,y,p1);
set(p1,'FaceColor','red','EdgeColor','none'); 
hold on
p2 = patch(isosurface(v,v,v,y,-S));
isonormals(v,v,v,y,p2);
set(p2,'FaceColor','blue','EdgeColor','none'); 
hold off
daspect([1 1 1]);
view(-30,30); 
camlight;
lighting phong
grid
axis([32 48 32 48 32 48])
set(gca,'fontsize',7)
title('3-D WAVELET ISOSURFACE (SEPARABLE TRANSFORM)')
set(gcf,'paperposition',[0.5 0.5 0 0]+[0 0 4 3])
print -djpeg95 dwt3D_plots
