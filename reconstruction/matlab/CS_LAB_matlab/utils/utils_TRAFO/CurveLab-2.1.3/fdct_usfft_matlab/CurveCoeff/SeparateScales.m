function S = SeparateScales(X,L)
% SeparatesScales: Separates an image into an (orthonormal) 
%                                    series of disjoint scales.
%  Usage:
%    S =  SeparateScales(X,L);
%  Inputs:
%    X     n by n image 
%    L     scale of the coarsest coefficients
%  Outputs: S data structure which contains the FT of X tapered at
%           dyadic scales.
%
% Description
%   For each s = 1, ..., J - L, 
%   S{s} is  a structure which contains the DFT of X tapered with a 
%   Meyer-style window which isolates frequencies near a dyadic 
%   frequency subband
%
%   S{s}         m by m matrix with m = 2*2^(j+1)); matrix
%                entries are Fourier coefficients of Img at 
%                -2^(j+1) <= k1,k2 < 2^(j+1) after windowing with a 
%                bandpass Meyer window which isolates frequencies 
%                in the range 2^(j-1) <= |k| <= 2^j.  
%                
%                j = L - 2 + s  
% See Also
%   Adj_ Separatescales, Curvelet02Xform
%
% By Emmanuel Candes, 2003-2004

  n = size(X,1);
  n2 = n/2;
  J = log2(n);
  deg = 3;
  
  scale = (L-1):(J-1);
  nscales = length(scale);   
  S = cell(1,nscales);  % Create data structure
    
  F = fft2_mid0(X)/sqrt(prod(size(X)));
  
  %      Partition at Coarse Level
  [ix,w] = CoarseMeyerWindow(L-1,deg);
  w = [0 reverse(w(2:length(w))) w]; 
  
  w2 = w'*w; % Build 2D-window
  
  l = 2^L;
  lx = (n2 - l + 1):(n2 + l); 
  S{1} = F(lx,lx).*w2;
  
  %      Loop to Get Partition for  j = L, ..., J - 1;
  for j = L:(J-1),
    
    l = min(2^(j+1),n2);
    wlo = zeros(1,l); whi = wlo;
	 
    [ixf, wf] = CoarseMeyerWindow(j-1,deg);
    wlo(ixf+1)  = wf;
    wlo = [0 reverse(wlo(2:length(wlo))) wlo]; 
    
    if j < J -1, 	    
      [ixp, wp] = CoarseMeyerWindow(j,deg); 
      whi(ixp+1) = wp; 
      whi = [0 reverse(whi(2:length(whi))) whi];
    end
    
    if j == J - 1, 
      whi = ones(1,2*l); 
    end
    
    w2 = sqrt((whi'*whi).^2 - (wlo'*wlo).^2); % Build 2D-window
    
    lx = (n2 - l + 1):(n2 + l); 
    S{j - L + 2} = F(lx,lx).*w2; 	  
  end
