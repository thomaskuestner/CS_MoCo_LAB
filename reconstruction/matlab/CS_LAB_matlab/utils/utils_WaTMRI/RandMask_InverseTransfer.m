function [mask] = RandMask_InverseTransfer(OMEGA, m, n)

M=zeros(m,n);
M(OMEGA)=1;
mask = ifftshift(M);
return