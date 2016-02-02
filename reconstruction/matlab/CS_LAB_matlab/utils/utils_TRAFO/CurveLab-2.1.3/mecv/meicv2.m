function x = meicv2(c,N1,N2,ns,nag)
% meicv2 - 2D inverse mirror-extended curvelet transform
% -----------------
% INPUT
% --
% c is a cell array which contains the curvelets coefficients. If
% tp=='ortho', then c{j}{l}(n1,n2) is the coefficient at scale j,
% direction l and spatial index (n1,n2). The directional index l
% iterates through the wedges in the first quadrant. Notice that, for
% the mirror-extended wave atoms, the spatial indices wrap around once.
% --
% N1, N2 are positive integers.
% --
% ns is the number of levels, including the coarsest level. ns =
% ceil(log2(min(N1,N2)) - 3) is commonly used.
% --
% nag is the number of angles used for the second coarsest level.
% nag is required to be a multiple of 4 and nag = 16 is often used.
% -----------------
% OUTPUT
% --
% x is an N1-by-N2 matrix. 
% -----------------
% Written by Lexing Ying and Laurent Demanet, 2007
  
  E1 = ceil(N1/3);  E2 = ceil(N2/3);  %E1 = 0;  E2 = 0;
  A1 = 2*(N1+E1);  A2 = 2*(N2+E2);
  fd = zeros(A1,A2);
  
  G1 = 4/3*N1;  G2 = 4/3*N2;
  
  for s=ns:-1:2
    R1 = 2^(s-ns)*G1;
    R2 = 2^(s-ns)*G2;
    
    idx1 = [ceil(-R1):floor(R1)];
    [wl,wr] = cvwindow((idx1+R1/1)/(R1/2));    tmpa = wl;
    [wl,wr] = cvwindow((idx1-R1/2)/(R1/2));    tmpb = wr;
    coef1 = tmpa.*tmpb;
    idx2 = [ceil(-R2):floor(R2)];
    [wl,wr] = cvwindow((idx2+R2/1)/(R2/2));    tmpa = wl;
    [wl,wr] = cvwindow((idx2-R2/2)/(R2/2));    tmpb = wr;
    coef2 = tmpa.*tmpb;
    lowpass = coef1'*coef2;
    
    idx1 = [ceil(-R1):floor(R1)];
    [wl,wr] = cvwindow((idx1+R1/2)/(R1/4));    tmpa = wl;
    [wl,wr] = cvwindow((idx1-R1/4)/(R1/4));    tmpb = wr;
    coef1 = tmpa.*tmpb;
    idx2 = [ceil(-R2):floor(R2)];
    [wl,wr] = cvwindow((idx2+R2/2)/(R2/4));    tmpa = wl;
    [wl,wr] = cvwindow((idx2-R2/4)/(R2/4));    tmpb = wr;
    coef2 = tmpa.*tmpb;
    tmppass = coef1'*coef2;
    hghpass = sqrt(1-tmppass.^2);
    
    pass = lowpass.*hghpass;
    [M1,M2] = size(pass);
    
    fh = zeros(M1,M2);
    %get fh
    
    nbangles = nag*2^(ceil((s-2)/2));
    %---------
    [M1,M2] = size(fh);
    nd = nbangles/4;
    cs = c{s};
    W1 = 2*R1/nd;  W2 = 2*R2/nd;

    %take only first quadrant
    cnt = 1;
    for g=nd/2:nd-1
      xs = R1/4-(W1/2)/4;    xe = R1;
      ys = -R2 + (2*g-1)*W2/2;		ye = -R2 + (2*g+3)*W2/2;
      xn = ceil(xe-xs);    yn = ceil(ye-ys);
      if(g==0)
        thts = atan2(-1.0, 1.0-1.0/nd);
        thtm = atan2(-1.0+1.0/nd, 1.0);
        thte = atan2(-1.0+3.0/nd, 1.0);
      elseif(g==nd-1)
        thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
        thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
        thte = atan2(1.0, 1.0-1.0/nd);
      else
        thts = atan2(-1.0+(2.0*g-1.0)/nd, 1.0);
        thtm = atan2(-1.0+(2.0*g+1.0)/nd, 1.0);
        thte = atan2(-1.0+(2.0*g+3.0)/nd, 1.0);
      end
      %fprintf(1,'%d %d %d\n',thts,thtm,thte);
      R21 = R2/R1;
      wpdata = fft2(cs{cnt}) / sqrt(numel(cs{cnt}));
      cnt = cnt+1;
      for xcur=ceil(xs):xe
        yfm = ceil( max([-R2, R21*xcur*tan(thts)]) );
        yto = floor( min([R2, R21*xcur*tan(thte)]) );
        ycur = yfm:yto;
        thtcur = atan2(ycur/R2,xcur/R1);
        [al,ar] = cvwindow((thtcur-thts)/(thtm-thts));
        [bl,br] = cvwindow((thtcur-thtm)/(thte-thtm));
        pou = al.*br;
        fh(mod(xcur,M1)+1,mod(ycur,M2)+1) = fh(mod(xcur,M1)+1,mod(ycur,M2)+1) + wpdata(mod(xcur,xn)+1,mod(ycur,yn)+1) .* pou;
      end
    end
    
    for f=nd-1:-1:nd/2
      ys = R2/4-(W2/2)/4;		  ye = R2;
      xs = -R1 + (2*f-1)*W1/2;		  xe = -R1 + (2*f+3)*W1/2;
      xn = ceil(xe-xs);		  yn = ceil(ye-ys);
      if(f==0)
        phis = atan2(-1.0, 1.0-1.0/nd);
        phim = atan2(-1.0+1.0/nd, 1.0);
        phie = atan2(-1.0+3.0/nd, 1.0);
      elseif(f==nd-1)
        phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
        phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
        phie = atan2(1.0, 1.0-1.0/nd);
      else
        phis = atan2(-1.0+(2.0*f-1.0)/nd, 1.0);
        phim = atan2(-1.0+(2.0*f+1.0)/nd, 1.0);
        phie = atan2(-1.0+(2.0*f+3.0)/nd, 1.0);
      end
      %fprintf(1,'%d %d %d\n',phis,phim,phie);
      R12 = R1/R2;
      wpdata = fft2(cs{cnt}) / sqrt(numel(cs{cnt}));
      cnt = cnt+1;
      for ycur=ceil(ys):ye
        xfm = ceil( max([-R1, R12*ycur*tan(phis)]) );
        xto = floor( min([R1, R12*ycur*tan(phie)]) );
        xcur = xfm:xto;
        phicur = atan2(xcur/R1, ycur/R2);
        [al,ar] = cvwindow((phicur-phis)/(phim-phis));
        [bl,br] = cvwindow((phicur-phim)/(phie-phim));
        pou = al.*br;
        fh(mod(xcur,M1)+1,mod(ycur,M2)+1) = fh(mod(xcur,M1)+1,mod(ycur,M2)+1) + wpdata(mod(xcur,xn)+1,mod(ycur,yn)+1) .* pou';
      end
    end
    
    %put back into fd
    fd(mod(idx1,A1)+1,mod(idx2,A2)+1) = fd(mod(idx1,A1)+1,mod(idx2,A2)+1) + pass .* fh(mod(idx1,M1)+1,mod(idx2,M2)+1);
  end
  
  if(1)
    s = 1;
    R1 = 2^(s-ns)*G1;
    R2 = 2^(s-ns)*G2;
    idx1 = [ceil(-R1):floor(R1)];
    [wl,wr] = cvwindow((idx1+R1/1)/(R1/2));    tmpa = wl;
    [wl,wr] = cvwindow((idx1-R1/2)/(R1/2));    tmpb = wr;
    coef1 = tmpa.*tmpb;
    idx2 = [ceil(-R2):floor(R2)];
    [wl,wr] = cvwindow((idx2+R2/1)/(R2/2));    tmpa = wl;
    [wl,wr] = cvwindow((idx2-R2/2)/(R2/2));    tmpb = wr;
    coef2 = tmpa.*tmpb;
    pass = coef1'*coef2;
    [M1,M2] = size(pass);
    
    cs = c{s};
    tmp = fft2(cs{1}) / sqrt(numel(cs{1}));
    tmp = tmp/4;
    tmp = mescatter(tmp,0);
    tmp = mescatter(tmp',0)';
    [K1,K2] = size(tmp);
    fh = zeros(M1,M2);
    fh(mod(idx1,M1)+1,mod(idx2,M2)+1) = tmp(mod(idx1,K1)+1,mod(idx2,K2)+1);
    
    fd(mod(idx1,A1)+1,mod(idx2,A2)+1) = fd(mod(idx1,A1)+1,mod(idx2,A2)+1) + pass .* fh(mod(idx1,M1)+1,mod(idx2,M2)+1);
  end
  
  fd = mecombine(fd,E1);
  fd = mecombine(fd',E2)';

  x = idct2(fd);
  
  