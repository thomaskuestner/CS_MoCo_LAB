function [K Kconv] = getK(G)
	K = G.KERNEL;
	if nargout > 1,
		Kconv = G.kernel;
	end
