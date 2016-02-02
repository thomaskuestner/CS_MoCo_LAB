function res = mtimes(A,x)

if A.adjoint == 0 %A*x
    res = dct(x);
    res = res(A.J);
else %At*x
    z = zeros(A.n,1);    
    z(A.J) = x;
    res = idct(z);
end
