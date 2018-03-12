function img = fdct_usfft_dispcoef(C)

% fdct_usfft_dispcoef - Returns an image containing all the curvelet coefficients
%
% Input
%     C         Curvelet coefficients 
%
% Output
%     img       Image containing all the curvelet coefficents. The coefficents are rescaled so that
%               the largest coefficent in each subband has unit norm.
%
  
  [m,n] = size(C{end}{1});
  nbscales = floor(log2(min(m,n)))-3;
  
  img = C{1}{1};  img = img/max(max(abs(img))); %normalize
  for sc=2:nbscales-1
    nd = length(C{sc})/4;
    wcnt = 0;
    
    ONE = [];
    for w=1:nd
      ONE = [ONE, C{sc}{wcnt+w}];
    end
    wcnt = wcnt+nd;
    
    TWO = [];
    for w=1:nd
      TWO = [TWO; C{sc}{wcnt+w}];
    end
    wcnt = wcnt+nd;
    
    THREE = [];
    for w=1:nd
      THREE = [C{sc}{wcnt+w}, THREE];
    end
    wcnt = wcnt+nd;
    
    FOUR = [];
    for w=1:nd
      FOUR = [C{sc}{wcnt+w}; FOUR];
    end
    wcnt = wcnt+nd;
    
    [p,q] = size(img);
    [a,b] = size(ONE);
    [g,h] = size(TWO);
    m = 2*a+g;    n = 2*h+b; %size of new image
    scale = max(max( max(max(abs(ONE))),max(max(abs(TWO))) ), max(max(max(abs(THREE))), max(max(abs(FOUR))) )); %scaling factor
    
    new = 0.5 * ones(m,n); %background value
    new(a+1:a+g,1:h) = FOUR /scale;
    new(a+g+1:2*a+g,h+1:h+b) = THREE /scale;
    new(a+1:a+g,h+b+1:2*h+b) = TWO /scale;
    new(1:a,h+1:h+b) = ONE /scale; %normalize
    
    dx = floor((g-p)/2);    dy = floor((b-q)/2);
    
    new(a+1+dx:a+p+dx,h+1+dy:h+q+dy) = img;
    
    img = new;
  end
