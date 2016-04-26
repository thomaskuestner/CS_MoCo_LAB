function y = cshift3D(x, m, d)

% 3D Circular Shift
%
% USAGE:
%   y = cshift3D(x, m, d)
% INPUT:
%   x - N1 by N2 by N3 array
%   m - amount of shift
%   d - dimension of shift (d = 1,2,3)
% OUTPUT:
%   y - array x will be shifed by m samples down
%       along dimension d
%
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

[N1, N2, N3] = size(x);
switch d
case 1
   n = 0:N1-1;
   n = mod(n-m, N1);
   y = x(n+1,:,:);
case 2
   n = 0:N2-1;
   n = mod(n-m, N2);
   y = x(:,n+1,:);
case 3
   n = 0:N3-1;
   n = mod(n-m, N3);
   y = x(:,:,n+1);
end

