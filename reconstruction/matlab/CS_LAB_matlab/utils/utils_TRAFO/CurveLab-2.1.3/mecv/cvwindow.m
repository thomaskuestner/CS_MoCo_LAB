function [wl,wr] = cvwindow(x)

  wr = zeros(size(x));
  wl = zeros(size(x));
  
  eps = 1e-16;
  sml = find(x<eps); %too small
  mid = find(x>=eps & x<=1-eps); %just right
  lrg = find(x>1-eps); %too large
  
  wl(sml) = 0;  wr(sml) = 1;
  wl(lrg) = 1;  wr(lrg) = 0;
  
  xmid = x(mid);
  a = exp(1-1./(1-exp(1-1./(1-xmid))));
  b = exp(1-1./(1-exp(1-1./xmid)));
  n = sqrt(a.^2 + b.^2);
  wl(mid) = a./n;  wr(mid) = b./n;
  
