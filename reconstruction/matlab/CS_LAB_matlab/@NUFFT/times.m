function ress = times(a,bb)
% performs only interpolation nufft

% for n=1:size(bb,3)
% b = bb(:,:,n);
% if a.adjoint
% 	b = b(:).*a.w(:);
% 	res = nufft_adj(b, a.st)/sqrt(prod(a.imSize));;
% 	res = reshape(res, a.imSize(1), a.imSize(2));
% 	res = fft2c(res);
% else
% 	b = reshape(b,a.imSize(1),a.imSize(2));
% 	b = ifft2c(b);
% 	res = (nufft(b, a.st)/sqrt(prod(a.imSize))).*a.w(:);
% 	res = reshape(res,a.dataSize(1),a.dataSize(2));
% end
% ress(:,:,n) = res;
% end
% 

ress = mtimes(a,bb);
