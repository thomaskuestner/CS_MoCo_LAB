function w = udwt(x, J, h0, h1)

% Undecimated Discrete Wavelet Transform
%
% INPUT
%   x : input signal
%   J : number of stages
%   h0 : low-pass analysis filter
%   h1 : high-pass analysis filter
% OUTPUT
%   w  : wavelet coefficients
%
% % Example:
%  [h0, h1, g0, g1] = daubf(3);
%  N = 20;
%  x = rand(1,N);
%  J = 3;
%  w = udwt(x, J, h0, h1);
%  y = iudwt(w,  J, g0, g1);
%  err = x - y(1:N);
%  max(abs(err))

R = sqrt(2);
h0 = h0/R;
h1 = h1/R;

w = cell(1,J);
for j = 1:J
    w{j} = upfirdn(h1, x, 2^(j-1), 1);
    x = upfirdn(h0, x, 2^(j-1), 1);
end
w{J+1} = x;

% Ivan Selesnick
% selesi@nyu.edu
% NYU - School of Engineering
