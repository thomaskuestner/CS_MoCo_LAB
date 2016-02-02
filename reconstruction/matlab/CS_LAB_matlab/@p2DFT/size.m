function [resx,resy] = size(a,n)

if a.adjoint
    res = a.imSize;
else
    res = a.dataSize;
end

if nargout<=1
        if exist('n')
                res = res(n);
        end
    resx = res;
    resy = [];
elseif nargout==2
    resx = res(1);
    resy = res(2);
else
   error('size return max of 2 values');
end

