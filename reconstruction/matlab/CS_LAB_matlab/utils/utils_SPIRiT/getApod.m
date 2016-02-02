function res = getApod(kernel,W,N,os)

% Given a kernel, grid size, kernel size and grid oversampling
% function will return the right apodization function

Gap = GRID(0,1,N,kernel,W,os);
res = crop(ifft2c(Gap.'*1),[N,N]);
