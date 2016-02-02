function x = Adj_USFFT(n,y,omega,D,L,center)
% USFFT_simple -- 1d adjoint unequispaced Fourier transform
% Usage:
%   x = USFFT(n,y,omega,D,L,center)
% Inputs:
%   n      length of output x
%   y	   vector of length m, 
%   omega  vector of sampled frequencies (length m)
%   center 0 for unbiased FT, 1 otherwise
% Outputs:
%   x      vector of length n
% Description:
%  Evaluates the adjoint of the FT
%   
%  If center = 0, 
%
%  x(t) = sum_{k} exp(i omega_k t) y(k), -n/2 <= t < n/2
%
%  If center = 1, 
%
%  x(t) = sum_{k} exp(i omega_k t) y(k), 0 <= t < n
%  
%  See Also
%   USFFT, USFT_simple, Evaluate_FT
%  Notes:
%   To work properly, D*n must be even
%
% By Emmanuel candes, 2003-2004

  if nargin < 6,
    center = 0;
  end
  
  if nargin < 5,
    L = 4;
  end
  
  if nargin < 4,
    D = 16;
  end
  
  n2 = n/2;
  N = D*n; N2 = N/2;
  m = length(omega);
  
  if rem(L,2) == 1,
    L = L + 1;
  end
  
  % Nearest point calculations
  
  near_index = round(N *omega/(2*pi));
  row = near_index + N/2 + 1;
  
  delta = omega - 2*pi*near_index/N; 
  
  % Make matrix y(k) * delta(k)^l, l = 0, 1, ..., L-1
  
  FA = repmat(y(:),1,L); 
  for l = 2:L;
    FA(:,l) = FA(:,l-1) .* delta/(l-1);
  end
  
  % Sum entries corresponding to the same index 
  
  F = zeros(N,L);
  for k = 1:m,
    F(row(k),:) = F(row(k),:) + FA(k,:);
  end
  
  % Take IFFT along columns
  X = ifft_mid0(F) * N; 
  X = X((N2-n2+1):(N2+n2),:);
  
  t = -n2:(n2-1);
  tp = ones(n,L); 
  
  for k = 2:L;
    tp(:,k) = tp(:,k-1) .* (i*t.');
  end	
  
  x = sum(X.*tp,2);
  
