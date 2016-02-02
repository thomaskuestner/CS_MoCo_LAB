% test_mecv2.m - this is a demo file which plots several mirror-extended curvelets at different scales

N = 512;
ns = ceil(log2(N) - 3);
nag = 16;

z = zeros(N,N);
c = mefcv2(z,N,N,ns,nag);

for is=ns-2:ns
  d = c;
  iw = floor(length(d{is})/2);
  [n1,n2] = size(d{is}{iw});
  d{is}{iw}(ceil(n1/4),1) = 1;
  y = meicv2(d,N,N,ns,nag);
  figure; imagesc(real(y)); axis xy; axis equal; axis tight; colormap(1-gray); colorbar;
  title(sprintf('curvelet near the boundary, j=%d, l=%d',is,iw));
  
  d = c;
  iw = floor(length(d{is})/2);
  [n1,n2] = size(d{is}{iw});
  d{is}{iw}(ceil(n1/4),ceil(n2/4)) = 1;
  y = meicv2(d,N,N,ns,nag);
  figure; imagesc(real(y)); axis xy; axis equal; axis tight; colormap(1-gray); colorbar;
  title(sprintf('curvelet at the center, j=%d, l=%d',is,iw));
end

