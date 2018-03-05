function a = subsasgn( a,s,b )
%SUBSASGN for nested cells
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

if(length(s) > 1)
    posCell = size(s(1).subs,2);
    posCell = s(1).subs{1,posCell};
    posMat = cell2mat(s(2).subs);
    [out,idx] = flattenCellMatrix(a.data);
    out{posCell}(sub2ind(size(out{posCell}),posMat(1),posMat(2),posMat(3))) = b;
    a.data = reconFlatCellMatrix(out,idx);
    
else
    if(strcmp(s.type,'{}'))
        pos = size(s.subs,2);
        pos = s.subs{1,pos};
        [out,idx] = flattenCellMatrix(a.data);
        out{pos} = b;
        a.data = reconFlatCellMatrix(out,idx);
        
    else
        pos = size(s.subs,2);
        pos = s.subs{1,pos};
        % a = a.data;
        % meta = a.meta;
        if(isa(b,'TRAFO'))
            if(size(s.subs,2) == 5)
                a.data(pos,:) = b.data; % for kernelImg (due to 5D dataset => cha-cha)
                if iscell(b.meta) ~= 1
                    a.meta = b.meta; %changed for curvelab!
                else
                    a.meta(pos,:) = b.meta;
                end
            else
                a.data(1,pos) = b.data;
                if(~iscell(b.meta) || length(b.meta) > 1)
                    a.meta = b.meta; %changed for curvelab!
                else
                    a.meta(1,pos) = b.meta;
                end
            end
        % a = TRAFO(a,meta);
        else
            [out,idx] = flattenCellMatrix(a.data(1,pos));
            for i=1:length(out)
                out{i}(:) = b;
            end
            a.data(1,pos) = reconFlatCellMatrix(out,idx);
        end
    end
end

end

