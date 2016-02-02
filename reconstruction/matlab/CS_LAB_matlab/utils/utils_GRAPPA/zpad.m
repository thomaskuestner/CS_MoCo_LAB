function res = zpad(x,sx,sy,sz,st)
% res = zpad(x,sx,sy,sz,st)
% Zero pad a matrix around its center
% or zero pad nested cell arrays
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

if(nargin < 2)
    error('zpad(): Must have a target size');
end

if(isa(x,'TRAFO'))
    if(nargin > 2)
        error('zpad(): Just two input arguments allowed for zeropadding of nested cell arrays');
    end
    if(~isa(sx,'TRAFO'))
        error('zpad(): Invalid input argument');
    end
    sxData = getData(sx);
    flatCell = flattenCellMatrix(getData(x));
    res = cellfun(@(x) zeros(x), sxData{2,1}, 'UniformOutput', false);
    for i=1:length(flatCell)
       res{1,i} = zpadKernel(flatCell{1,i}, sxData{2,1}{1,i}); 
    end
    res = reconFlatCellMatrix(res, sxData{1,1});
    res = TRAFO(res,getMeta(sx));
else
    if(nargin == 2)
        s = sx;
    elseif(nargin == 3)
        s = [sx,sy];
    elseif(nargin == 4)
        s = [sx,sy,sz];
    elseif(nargin == 5)
        s = [sx,sy,sz,st];
    end
    
    res = zpadKernel(x,s);
end

    function res = zpadKernel(x,s)
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
        % modifications by Thomas Kuestner 2013


        m = size(x);
        if length(m) < length(s)
            m = [m, ones(1,length(s)-length(m))];
        end

        if sum(m==s)==length(m)
            res = x;
            return;
        end

        res = zeros(s);

        idx = cell(1,length(s));
        for n=1:length(s)
            if(mod(s(n),2) == 0)
                idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
            else
                idx{n} = floor(s(n)/2)+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2)-1;
            end
            helper = [idx{n}(1) <= 0, idx{n}(end) > s(n)];
            if(any(helper))
                if(all(helper)), error('zpad(): Both index out of bounds'); end;
                hShift = [abs(idx{n}(1)) + 1, idx{n}(end) - s(n)];
                op = {idx{n}(1), idx{n}(end); '+', '-'; '<', '>'; s(n) + 1, 0};
                eval(sprintf('if(op{1,~helper} %s hShift(helper) %s %d), idx{n} = idx{n} %s hShift(helper); else idx{n} = idx{n}(idx{n} %s %d);end;',op{2,helper}, op{3,helper}, op{4,helper}, op{2,helper}, op{3,~helper}, op{4,~helper}));
           end
        end

        % this is a dirty ugly trick
        cmd = 'res(idx{1}';
        for n=2:length(s)
            cmd = sprintf('%s,idx{%d}',cmd,n);
        end
        cmd = sprintf('%s)=x;',cmd);
        eval(cmd);
    
    end

end





