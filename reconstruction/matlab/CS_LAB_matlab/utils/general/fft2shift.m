function res = fft2shift(x)

warning('Deprecated! Consider using fftnshift(x) instead');
res = fftnshift(x,1:2);

end


