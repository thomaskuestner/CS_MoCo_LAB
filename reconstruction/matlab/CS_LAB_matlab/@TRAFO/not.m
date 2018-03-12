function out = not( in )
%NOT for nested cells
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

[out, idx] = flattenCellMatrix(in.data);
out = cellfun(@(x) ~x, out, 'UniformOutput', false);
out = reconFlatCellMatrix(out,idx);
out = TRAFO(out,in.meta);

end

