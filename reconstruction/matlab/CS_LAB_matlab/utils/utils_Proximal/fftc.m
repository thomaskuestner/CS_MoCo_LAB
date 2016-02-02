function res = fftc(x,dim)
% res = fftc(x,dim)
res = 1/sqrt(size(x,dim))*fftshift(fft(ifftshift(x,dim),[],dim),dim);

