function [h0, h1, g0, g1] = daubf(K,str);
% [h0, h1, g0, g1] = daubf(K);
% h0 - low-pass analysis
% h1 - high-pass analysis
% g0 - low-pass analysis
% g1 - high-pass analysis
%
% K zeros at z=-1
% Use daubf(K,'mid') for mid-phase type

% Ivan Selesnick
% selesi@nyu.edu
% NYU - School of Engineering

[h,s,g] = maxflatI(K,K-1);

r = roots(g);
r = r(abs(r) < 1);
q = real(poly(r));

if nargin > 1
	if strcmp(str,'mid')
		q = sfact_mid(g);
	end
end

q = q/sum(q);                 % normalize
h0 = q;                       % set  h0 = q;
for k = 1:K                   % make h0 = q * [(z^(-1)+1)/2]^K
   h0 = conv(h0,[1 1]/2);
end
h0 = sqrt(2)*h0;              % normalize so that sum(h0) = sqrt(2)

h1 = h0(end:-1:1);
h1(2:2:end) = -h1(2:2:end);

g0 = h0(end:-1:1);
g1 = h1(end:-1:1);
