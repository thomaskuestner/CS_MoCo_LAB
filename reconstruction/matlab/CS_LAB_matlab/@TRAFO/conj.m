function out = conj( in )
%CONJ for nested cells
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

[out, idx] = flattenCellMatrix(in.data);
out = cellfun(@(x) conj(x), out, 'UniformOutput', false);
out = reconFlatCellMatrix(out,idx);
out = TRAFO(out,in.meta);

end

