 function tf = streq(a,b,n)
%function tf = streq(a, b [,n])
%	return 1 if two strings are equal (optionally only up to 1st n chars)
if nargin == 2
	tf = strcmp(a,b);
elseif nargin == 3
	tf = strncmp(a,b,n);
else
	help(mfilename), error(mfilename)
end
