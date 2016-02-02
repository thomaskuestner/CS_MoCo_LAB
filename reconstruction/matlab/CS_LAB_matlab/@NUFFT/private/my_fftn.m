function X = my_fftn(x)

X = fft(x);
dd = ndims(x);
for id=2:dd
	ord = 1:dd;
	ord([1 id]) = [id 1];
	x = permute(X, ord);
	X = fft(x);
	X = ipermute(X, ord);
end
