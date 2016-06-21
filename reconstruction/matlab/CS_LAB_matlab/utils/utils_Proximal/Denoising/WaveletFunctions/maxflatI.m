function [h,h1,h2] = maxflatI(K,M)
% [h,h1,h2] = maxflatI(K,M)
% Maxflat Type-I FIR filter 
%   2K zeros at z=-1
%   2M zeros away from z=-1
%   h = conv(h1,h2); 
%   h1 : all zeros at z=-1
%   h2 : all other zeros
%
% Note: if K = M+1, then h is halfband.
%
% Reference:
% O. Herrmann, "On the approximation problem in Nonrecursive
% Digital Filter Design", IEEE Trans. on Circuit Theory,
% Vol. 18, No. 3, May 1971, pp. 411-413
%
% % Example
% [h,h1,h2] = maxflatI(4,6);

% Ivan Selesnick
% selesi@nyu.edu
% NYU - School of Engineering


h2 = 1;
hi = 1;
c  = 1;
for k = 1:M
   hi = conv(hi,[-1 2 -1]/4);
   c  = c*(K-1+k)/k;
   h2 = [0 h2 0] + c*hi;
end

h1 = 1;
for k = 1:2*K
   h1 = conv(h1,[1 1]/2);
end

h = conv(h1,h2);

