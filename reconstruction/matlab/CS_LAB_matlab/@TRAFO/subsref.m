function out = subsref( in, s )
%SUBSREF for nested cells
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

% return array with concatenated data values
if(strcmp(s.type,'.'))
    out = in;
else    
    if(size(s.subs,2) == 1)
        out = flattenCellMatrix(in.data);
        if(isnumeric(s.subs{1}))
            out = cellfun(@(x) x(s.subs{:}), out.', 'UniformOutput', false); % for totEnergy(c)
            out = TRAFO(out, in.meta);
        else
            out = cell2mat(cellfun(@(x) x(s.subs{:}), out.', 'UniformOutput', false)); % for helper(:)
        end
    else
        pos = size(s.subs,2);
        out = in.data;
        if(pos == 5)
            if iscell(in.meta) ~= 1
                out = TRAFO(out(s.subs{1,pos},:),in.meta); % for kernelImg(:,:,:,:,i) (curvelab)
            else
                out = TRAFO(out(s.subs{1,pos},:),in.meta(s.subs{1,pos},:)); % for kernelImg(:,:,:,:,i)
            end
        else
            if iscell(in.meta) ~= 1
                out = TRAFO(out(1,s.subs{1,pos}),in.meta); % for 4D datasets (curvelab)
            else 
                out = TRAFO(out(1,s.subs{1,pos}),in.meta(1,s.subs{1,pos})); % for 4D datasets
            end
        end
    end
end
end

