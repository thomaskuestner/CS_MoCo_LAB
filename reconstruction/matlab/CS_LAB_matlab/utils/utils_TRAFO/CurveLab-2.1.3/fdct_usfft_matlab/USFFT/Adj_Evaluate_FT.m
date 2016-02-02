function x = Adj_Evaluate_FT(y,shift,boxlen,center,w)
% Adj_Evaluate_FT -- 1d unequispaced adjoint Fourier transform
% (with windowing). Adjoint of Evaluate_FT.  
% Usage:
%   x = Evaluate_FT(y,shift,boxlen,center,w)
% Inputs:
%   y	      vector of length 2n
%   shift     vector of shifts
%   boxlen    half length of interval around each value of shift
%   center    boolean variable: 0/1 = unbiased/biased FT
%   w         tapering window of size 2*boxlen
% Outputs:
%   x      vector of length: n
% Description: 
%  For w = 1, evaluates the Adj FT of a signal irregularly
%  sampled in frequency omega(j,k) = shift(j) + k, -boxlen <= k
%  <boxlen If center = 0,
%
%  x(t) = sum_m exp(i 2pi omega(m) t/n) y(k), -n/2 <= t < n/2
%
%  If center = 1, 
%
%  x(t) = sum_m exp(i 2pi omega(m) t/n) y(k), 0 <= t < n
% 
%  For arbitrary windows, y is first tapered with w.
%
%  Adj_Evaluate_FT uses two different methods for computing the
%  unesquispaced Fourier transform. If the number of shifts is less
%  than 8, it uses an exact method, otherwise it uses the
%  accelerated Adjoint USFFT. 
%  See Also
%    Evaluate_FT, Adj_USFFT, Adj_USFT_simple
%
% By Emmanuel candes, 2003-2004

if nargin < 5,
  w = ones(1,2*boxlen);
end

if nargin < 4,
  center = 0;
end

n  = length(y)/2;
boxcnt = length(shift);

if (boxcnt <= 16)
  x = Adj_USFT_simple(y,shift,boxlen,center,w);
else
  % Make frequency grid
  k = (-boxlen):(boxlen -1);
  omega = meshgrid(shift,k) + meshgrid(k,shift).';
  omega = 2*pi*(omega(:))./n;
  
  y = repmat(w(:),boxcnt,1) .* y; % Tapering
  x = Adj_USFFT(n,y,omega,16,6,center)./sqrt(n); % Take Adjoint USFFT
end
