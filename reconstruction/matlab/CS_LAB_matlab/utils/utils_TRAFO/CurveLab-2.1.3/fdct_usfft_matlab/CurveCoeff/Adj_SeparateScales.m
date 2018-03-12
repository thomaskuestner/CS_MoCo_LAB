function X = Adj_SeparateScales(S,L)
% SeparatesScales: Adjoint of the scale separation operation (also
%                          its inverse)
%  Usage:
%    X =  Adj_SeparateScales(S,L);
%  Inputs:
%    S data structure which contains the FT of X at dyadic scales. 
%    L     scale of the coarsest coeffcients
%  Outputs:
%    X n by n image
%  Description
%    Multiplies frequency sample in each dyadic corona by a Meyer window 
%    and adds the contributions from all the different scales. 
% See Also
%   Separatescales, Adj_Curvelet02Xform
%
% By Emmanuel Candes, 2003-2004

  number_scales = length(S);   
  n = size(S{number_scales},1);
  n2 = n/2;
  J = log2(n);
  deg = 3;
  F = zeros(n);
  
  %      Unfold Partition at Coarse Level
  
  [ix,w] = CoarseMeyerWindow(L-1,deg);
  w = [0 reverse(w(2:length(w))) w];
  
  w2 = w'*w; % Build 2D-window
  
  l = 2^L;
  lx = (n2 - l + 1):(n2 + l); % ly = lx;
  F(lx,lx) = S{1}.*w2;
  
  %      Loop to Unfold Partition for  j = L, ..., J - 1;
  
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
    
    lx = (n2 - l + 1):(n2 + l); % ly = lx;
    F(lx,lx) = F(lx,lx) + S{j-L+2}.*w2;
  end
  
  X  = ifft2_mid0(F)*sqrt(prod(size(F)));
