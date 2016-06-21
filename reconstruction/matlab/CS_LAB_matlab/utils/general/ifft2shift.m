function res = ifft2shift(x)

warning('Deprecated! Consider using ifftnshift(x) instead');
res = ifftnshift(x,1:2);

end
