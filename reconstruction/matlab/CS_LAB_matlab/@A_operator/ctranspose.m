function At = ctranspose(A)
At.func = A.transpose;
At.transpose = A.func;
At = class(At,'A_operator');