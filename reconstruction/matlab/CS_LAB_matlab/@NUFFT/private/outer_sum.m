 function ss = outer_sum(xx,yy)
%function ss = outer_sum(xx,yy)
%	compute an "outer sum" x + y'
%	that is analogous to the "outer product" x * y'
%	Input
%		xx [nx,1]
%		yy [1,ny]
%	Output
%		ss [nx,ny]	ss(i,j) = xx(i) + yy(j)
if ~nargin, help outer_sum, end

	xx = xx(:);
	yy = yy(:)';
	nx = length(xx);
	ny = length(yy);
	xx = xx(:, ones(1,ny));
	yy = yy(ones(1,nx),:);
%	[xx, yy] = ndgrid(xx, yy);
	ss = xx + yy;
