function  res = ctranspose(A)

A.adjoint = xor(A.adjoint,1);
res = A;