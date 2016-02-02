function ress = mtimes(a,bb)
% performs the normal nufft

for n=1:size(bb,3)
b = bb(:,:,n);

if a.adjoint
	b = b(:).*a.w(:);
	res = nufft_adj(b, a.st)/sqrt(prod(a.imSize));
	res = reshape(res, a.imSize(1), a.imSize(2));

else
	b = reshape(b,a.imSize(1),a.imSize(2));
	res = nufft(b, a.st)/sqrt(prod(a.imSize)).*a.w(:);
	res = reshape(res,a.dataSize(1),a.dataSize(2));
end
ress(:,:,n) = res;
end

