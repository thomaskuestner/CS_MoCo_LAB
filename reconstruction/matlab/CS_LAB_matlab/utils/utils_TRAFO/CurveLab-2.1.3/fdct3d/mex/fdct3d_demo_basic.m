disp(' ');
disp(['fdct3d_demo_basic.m -- This demo displays a 3d curvelet'])
disp (['both in the spatial and frequency domains.']);
disp(' ');
disp(['This is achieved by setting all the coefficients in the curvelet'])
disp(['domain to zero except that at the required location (which'])
disp(['is set to one). The curvelet is obtained by taking the'])
disp(['adjoint curvelet transform. Notice how the curvelet is sharply '])
disp(['localized in both space and frequency.']); 
disp(' ');
disp('The curvelet displayed is at the 4th scale. Its energy in the frequency')
disp('domain is concentrated in the direction of (x,y,z)=(1,-1,-1). In the spatial')
disp('domain, it looks like a disc with the normal direction equal to (1,-1,-1).')
disp('The curvelet oscillates in the normal direction.')
disp(' ');

% fdct3d_demo_basic.m -- This demo displays a curvelet both in the spatial and frequency domains.

m = 128;  n = 128;  p = 128;
s = 4;
w = 1;

X = zeros(m,n,p);

disp('Take forward 3d curvelet transform:');
tic; C = fdct3d_forward(X); toc; 

[t1,t2,t3] = size(C{s}{w});
t1 = ceil((t1+1)/2);  t2 = ceil((t2+1)/2);  t3 = ceil((t3+1)/2);
C{s}{w}(t1,t2,t3) = 1;

disp('Take inverse 3d curvelet transform:');
tic; Y = fdct3d_inverse(C); toc; 

F = ifftshift(fftn(Y));
%reorder the data for display
Y = permute(Y, [2,1,3]);
F = permute(F, [2,1,3]);

%display 1
h = slice(real(Y),m/2,n/2,p/2);
alpha('color')
set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
alphamap('rampdown')
alphamap('increase',.1)
colormap(hsv)
set(gcf,'Renderer','zbuffer'); lighting phong
xlabel('x'); ylabel('y'); zlabel('z');
axis([0,m,0,n,0,p]);  colormap gray; axis equal;

disp(' ');
disp('Sliced display of the curvelet');
disp('Press any key to start the animation which displays the curvelet at different slices'); 
pause;

%display 2
for ix=m-1:-4:1
  subplot(1,2,1);
  h = slice(real(Y),ix,[],[]);
  alpha('color')
  set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
  alphamap('rampdown')
  alphamap('increase',.1)
  colormap(hsv)
  set(gcf,'Renderer','zbuffer'); lighting phong
  xlabel('x'); ylabel('y'); zlabel('z');
  axis([0,m,0,n,0,p]);  colormap gray;
  title('a curvelet: spatial viewpoint');
  
  subplot(1,2,2);
  h = slice(abs(F),ix,[],[]);
  alpha('color')
  set(h,'EdgeColor','none','FaceColor','interp','FaceAlpha','interp')
  alphamap('rampdown')
  alphamap('increase',.1)
  colormap(hsv)
  set(gcf,'Renderer','zbuffer'); lighting phong
  xlabel('x'); ylabel('y'); zlabel('z');
  axis([0,m,0,n,0,p]);  colormap gray;
  title('a curvelet: frequency viewpoint');
  
  pause(0.05);
end
