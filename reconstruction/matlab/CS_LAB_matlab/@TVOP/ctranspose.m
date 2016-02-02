function res = ctranspose(a)
a.adjoint = xor(a.adjoint,1);
res = a;

