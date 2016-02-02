function res = crop(x,sx,sy,sz,st)
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

if(nargin < 2)
	error('crop(): Must have a target size');
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
       res{1,i} = cropKernel(flatCell{1,i}, sxData{2,1}{1,i}); 
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
    
    res = cropKernel(x,s);
end


    function res = cropKernel(x,s)
        %  res = crop(x,sx,sy)
        %  crops a 2D matrix around its center.
        %
        %
        %  res = crop(x,sx,sy,sz,st)
        %  crops a 4D matrix around its center
        %
        %  
        %  res = crop(x,[sx,sy,sz,st])
        %  same as the previous example
        %
        %
        %
        %
        % (c) Michael Lustig 2007
        % modifications by Thomas Kuestner 2013


        m = size(x);
        if length(s) < length(m)
            s = [s, ones(1,length(m)-length(s))];
        elseif(length(s) > length(m))
            if(s(end) == 1)
                s = s(1:end-1);
            end
        end

        if sum(m==s)==length(m)
            res = x;
            return;
        end

        idx = cell(1,length(s));
        for n=1:length(s)
            if(mod(s(n),2) == 0)
                idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
            else
                idx{n} = floor(m(n)/2)+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2)-1;
            end
            
            helper = [idx{n}(1) <= 0, idx{n}(end) > m(n)];
            if(any(helper))
                if(all(helper)), error('crop(): Both index out of bounds'); end;
                hShift = [abs(idx{n}(1)) + 1, idx{n}(end) - m(n)];
                op = {idx{n}(1), idx{n}(end); '+', '-'; '<', '>'; m(n) + 1, 0};
                eval(sprintf('if(op{1,~helper} %s hShift(helper) %s %d), idx{n} = idx{n} %s hShift(helper); else idx{n} = idx{n}(idx{n} %s %d);end;',op{2,helper}, op{3,helper}, op{4,helper}, op{2,helper}, op{3,~helper}, op{4,~helper}));
                % for left boundary: 
                % if(op{1,2} + hShift(1) < m(n) + 1)
                %   idx{n} = idx{n} + hShift(1);
                % else
                %   idx{n} = idx{n}(idx{n} > 0);
                % end
                % for right boundary:
                % if(op{1,1} - hShift(2) > 0)
                %   idx{n} = idx{n} - hShift(2);
                % else
                %   idx{n} = idx{n}(idx{n} < m(n) + 1);
                % end
           end
        end

        % this is a dirty ugly trick
        cmd = 'res = x(idx{1}';
        for n=2:length(s)
            cmd = sprintf('%s,idx{%d}',cmd,n);
        end
        cmd = sprintf('%s);',cmd);
        eval(cmd);
    end
    
end





