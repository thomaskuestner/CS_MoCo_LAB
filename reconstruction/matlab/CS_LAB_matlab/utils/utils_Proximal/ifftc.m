function res = ifftc(x,dim)
%res = ifftc(x,dim)
res = sqrt(size(x,dim))*fftshift(ifft(ifftshift(x,dim),[],dim),dim);

