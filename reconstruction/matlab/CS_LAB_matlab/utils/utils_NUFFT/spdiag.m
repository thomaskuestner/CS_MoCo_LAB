 function b = spdiag(a)
%function b = spdiag(a)
%	create a sparse matrix with diagonal given by a

	a = a(:);
	b = spdiags(a, 0, length(a), length(a));
