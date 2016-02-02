function res = zpad(x,sx,sy,sz,st)
%  res = zpad(x,sx,sy)
%  Zero pads a 2D matrix around its center.
%
%
%  res = zpad(x,sx,sy,sz,st)
%  Zero pads a 4D matrix around its center
%
%  
%  res = zpad(x,[sx,sy,sz,st])
%  same as the previous example
%
%
% (c) Michael Lustig 2007

if nargin < 2
	error('must have a target size')
end

if nargin == 2
	s = sx;
end

if nargin == 3
    s = [sx,sy];
end

if nargin == 4
    s = [sx,sy,sz];
end

if nargin == 5
    s = [sx,sy,sz,st];
end

    m = size(x);
    if length(m) < length(s)
	    m = [m, ones(1,length(s)-length(m))];
    end
	
    if sum(m==s)==length(m)
	res = x;
	return;
    end

    res = zeros(s);
    
    for n=1:length(s)
	    idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
    end

    % this is a dirty ugly trick
    cmd = 'res(idx{1}';
    for n=2:length(s)
    	cmd = sprintf('%s,idx{%d}',cmd,n);
    end
    cmd = sprintf('%s)=x;',cmd);
    eval(cmd);





