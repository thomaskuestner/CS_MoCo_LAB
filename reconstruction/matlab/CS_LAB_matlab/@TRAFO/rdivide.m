function out = rdivide( dividend, divisor )
%RDIVIDE for nested cells
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

% determine input types
inType = {class(dividend), class(divisor)};
if(strcmp(inType{1},'TRAFO') && strcmp(inType{2},'TRAFO'))
    [outA, idx] = flattenCellMatrix(dividend.data);
    [outB, idxB] = flattenCellMatrix(divisor.data); 
    meta = dividend.meta;
    meta_b = divisor.meta;
    clear 'dividend' 'divisor'
    
    % idx must be consistent
    if(~isempty(idx) && ~isempty(idxB) && ~all(cellfun(@(x,y) isequal(x,y), idx, idxB)))
        error('TRAFO::rdivide: Nested cells must have the same size');
    end
    if(iscell(meta) && iscell(meta_b))
        if(~isempty(meta) && ~isempty(meta_b))
            if(iscell(meta{1}) && iscell(meta_b{1}))
                if(~all(cellfun(@(x,y) isequal(x,y), flattenCellMatrix(meta), flattenCellMatrix(meta_b))))
                    error('TRAFO::rdivide: Unequal reconstruction information');
                end
            else
                if(~all(cellfun(@(x,y) isequal(x,y), meta, meta_b)))
                    error('TRAFO::rdivide: Unequal reconstruction information');
                end
            end
        end
    else 
        if ~isequal(meta, meta_b)
            error('TRAFO::rdivide: Unequal reconstruction information');
        end
    end
    out = cellfun(@(x,y) x./y, outA, outB, 'UniformOutput', false);

elseif(strcmp(inType{1},'TRAFO')) % assume divisor is double/uint/int
    [out, idx] = flattenCellMatrix(dividend.data);
    meta = dividend.meta;
    clear 'dividend'
    
    out = cellfun(@(x) x./divisor, out, 'UniformOutput', false);
    
elseif(strcmp(inType{2},'TRAFO'))
    [out, idx] = flattenCellMatrix(divisor.data);
    meta = divisor.meta;
    clear 'divisor'
    
    out = cellfun(@(x) dividend./x, out, 'UniformOutput', false);
    
else
    error('TRAFO::rdivide: Impossible constellation!');
end

out = reconFlatCellMatrix(out,idx);
out = TRAFO(out,meta); % return TRAFO object again without modifying input

end

