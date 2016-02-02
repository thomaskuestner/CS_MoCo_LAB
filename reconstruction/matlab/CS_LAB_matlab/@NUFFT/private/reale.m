 function y = reale(x, arg2)
%function y = reale(x)
%function y = reale(x, tol)
%function y = reale(x, 'warn')
%function y = reale(x, 'error')
%function y = reale(x, 'report')
%function y = reale(x, 'disp')
% return real part of complex data (with error checking).
% checks that imaginary part is negligible (or warning etc. if not)
% Copyright Jeff Fessler, The University of Michigan

com = 'error';
tol = 1e-13;
if nargin > 1
	if ischar(arg2)
		com = arg2;
	elseif isnumeric(arg2)
		tol = arg2;
	end
end

if streq(com, 'disp')
	;
elseif streq(com, 'warn')
	onlywarn = 1;
elseif streq(com, 'error')
	onlywarn = 0;
elseif streq(com, 'report')
	;
else
	error(sprintf('bad argument %s', com))
end

if max(abs(x(:))) == 0
	if any(imag(x(:)) ~= 0)
		error 'max real 0, but imaginary!'
	else
		y = real(x);
		return
end
end

frac = max(abs(imag(x(:)))) / max(abs(x(:)));
if streq(com, 'report')
	printf('imaginary part %g%%', frac * 100)
	return
end

if frac > tol
	t = sprintf('%s: imaginary fraction %g', mfilename, frac);
	if streq(com, 'disp')
		disp(t)
	elseif onlywarn
		disp(t)
	else
		error(t)
	end
end
y = real(x);
