% cplxdual3D_plots
% DISPLAY 3D WAVELETS OF CPLXDUAL3D.M

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
J = 4;
L = 2*2^(J+1);
N = L/2^J;
x = zeros(L,L,L);
w = cplxdual3D(x, J, Faf, af);
w{J}{1}{1}{1}{3}(N/2,N/2,N/2) = 1;
y = icplxdual3D(w, J, Fsf, sf);
figure(1)
clf
v = 1:L;
S = 0.0010;
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
axis([12 38 12 38 12 38])
set(gca,'fontsize',7)
title('3-D WAVELET ISOSURFACE (REAL PART)')
set(gcf,'paperposition',[0.5 0.5 0 0]+[0 0 3 3])
print -djpeg95 cplxdual3D_plots_1
print -depsc cplxdual3D_plots_1
yr = y;

w{J}{1}{1}{1}{3}(N/2,N/2,N/2) = 0;
w{J}{2}{2}{2}{3}(N/2,N/2,N/2) = 1;
y = icplxdual3D(w, J, Fsf, sf);
figure(2)
clf
v = 1:L;
S = 0.0010;
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
axis([12 38 12 38 12 38])
set(gca,'fontsize',7)
title('3-D WAVELET ISOSURFACE (IMAGINARY PART)')
set(gcf,'paperposition',[0.5 0.5 0 0]+[0 0 3 3])
print -djpeg95 cplxdual3D_plots_2
print -depsc cplxdual3D_plots_2
yi = y;

y = sqrt(yr.^2 + yi.^2);
figure(3)
clf
v = 1:L;
S = 0.0010;
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
axis([12 38 12 38 12 38])
set(gca,'fontsize',7)
title('3-D WAVELET ISOSURFACE (MAGNITUDE)')
set(gcf,'paperposition',[0.5 0.5 0 0]+[0 0 3 3])
print -djpeg95 cplxdual3D_plots_3
print -depsc cplxdual3D_plots_3

