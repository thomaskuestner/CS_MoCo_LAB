function res = fftc(x,dim)
% res = fftc(x,dim)
if nargin < 2
    dim = 2;
end

res = 1/sqrt(size(x,dim))*fftshift(fft(ifftshift(x,dim),[],dim),dim);

