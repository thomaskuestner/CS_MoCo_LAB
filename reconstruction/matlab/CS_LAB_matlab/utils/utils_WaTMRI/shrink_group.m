function z = shrink_group(x,y,nonneg,groups,g)
%
%  Group-wise or row-wise shrinkage operator: Shrink(x,y)
%

if nonneg  % nonnegativity
    x = max(x,0);
end

if nargin == 5  % group-wise
    t = max(0,1-y./sqrt(g*x.^2));
    z = t(groups).*x;
end

if nargin == 3  % row-wise
    t = max(0,1-y./sqrt(sum(x.^2,2)));
    z = bsxfun(@times,t,x);
end
end