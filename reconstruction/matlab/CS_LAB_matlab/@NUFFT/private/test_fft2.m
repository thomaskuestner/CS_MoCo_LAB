% test_fft2.m
% compare fft2 vs fft(fft())
%
% results on G5:
% fftn 0.147631
% fft2 0.133448
% fft1 0.130262
% fft0 0.0802386
%
% results on ir7:
% fftn 0.184001
% fft2 0.184357
% fft1 0.192388
% fft0 0.177634

N = 2^9;
M = 5;
rand('state', 0)
x = rand(N);

tic
for m=1:M, X1 = fftn(x); end
printf('fftn %g', toc/M)

tic
for m=1:M, Xm = my_fftn(x); end
printf('fftm %g', toc/M)

tic
for m=1:M, X4 = fft(fft(x).').'; end
printf('fft0 %g', toc/M)

tic
for m=1:M, X2 = fft2(x); end
printf('fft2 %g', toc/M)

tic
for m=1:M, X3 = fft(fft(x,[],2),[],1); end
printf('fft1 %g', toc/M)

[max_percent_diff(X1,X2) max_percent_diff(X1,X3) ...
max_percent_diff(X1,X4) max_percent_diff(X1,Xm)]

