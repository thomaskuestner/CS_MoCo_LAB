function y = iudwt(w, J, g0, g1)

% y = iudwt(w, J, g0, g1)
%
% Inverse undecimated discrete wavelet transform
% INPUT
%   w : wavelet coefficients
%   J : number of stages
%   g0 : low-pass synthesis filter
%   g1 : high-pass synthesis filter
% OUTPUT
%   y : output signal

R = sqrt(2);
g0 = g0/R;
g1 = g1/R;
N = length(g0) + length(g1);

y = w{J+1};
for j = J:-1:1
    M = 2^(j-1);
    lo = upfirdn(g0, y, M, 1);
    hi = upfirdn(g1, w{j}, M, 1);

    % Add signals, remove leading/trailing zeros
    Nhi = length(hi);
    L = M*(N/2-1);
    y = lo(L+1:Nhi-L) + hi(L+1:Nhi-L);
end

% Ivan Selesnick
% selesi@nyu.edu
% NYU - School of Engineering


