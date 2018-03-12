function out = ne( a,b )
%ne (~=) for nested cell arrays
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

% determine input types
inType = {class(a), class(b)};
if(strcmp(inType{1},'TRAFO') && strcmp(inType{2},'TRAFO'))  
    [outA, idx] = flattenCellMatrix(a.data);
    [outB, idxB] = flattenCellMatrix(b.data);
    meta = a.meta;
    meta_b = b.meta;
    clear 'a' 'b'
    out = false;

    if(~all(cellfun(@(x,y) isequal(x,y), idx, idxB)))
        error('TRAFO::ne: Nested cells must have the same size');
    end
    if (iscell(meta) && iscell(meta_b))
        if(~isempty(meta) && ~isempty(meta_b)) % unequal reconstruction information can happen -> this is what we are interested in to check
            if(iscell(meta{1}) && iscell(meta_b{1}))
                if(~all(cellfun(@(x,y) isequal(x,y), flattenCellMatrix(meta), flattenCellMatrix(meta_b))))
                    out = true;
    %                 error('TRAFO::ne: Unequal reconstruction information');
                end
            else
                if(~all(cellfun(@(x,y) isequal(x,y), meta, meta_b)))
                    out = true;
    %                 error('TRAFO::ne: Unequal reconstruction information');
                end
            end
        end
        out = all(cellfun(@(x,y) ~isequal(x,y), outA, outB)) | out;
    else
        if(~isequal(meta, meta_b))
            out = true; % can not be the same, because the size of the original image stored as meta info is different
        else
            out = all(cellfun(@(x,y) ~isequal(x,y), outA, outB)) | out; %check if entries are the same
        end
    end

elseif(strcmp(inType{1},'TRAFO')) % assume b is double/uint/int
    [out, ~] = flattenCellMatrix(a.data);
%     meta = a.meta;    
    clear 'a'
    
    out = all(cellfun(@(x) ~isequal(x,b), out));
    
elseif(strcmp(inType{2},'TRAFO'))
    [out, ~] = flattenCellMatrix(b.data);
%     meta = b.meta;
    clear 'b'
    
    out = all(cellfun(@(x) ~isequal(a,x), out));
    
else
    error('TRAFO::ne: Impossible constellation!');
end

end

