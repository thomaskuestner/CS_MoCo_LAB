function  res = partialDCT(n,m,J)

res.adjoint = 0;
res.n = n;
res.m = m;
res.J = J(1:m);

% Register this variable as a partialDCT class
res = class(res,'partialDCT');
