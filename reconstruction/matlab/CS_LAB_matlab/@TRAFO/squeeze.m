function out = squeeze( in )
%SQUEEZE nested cell arrays
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

[out, idx] = flattenCellMatrix(in.data);
out = cellfun(@(x) squeeze(x), out, 'UniformOutput', false);
out = reconFlatCellMatrix(out, idx);
out = TRAFO(out,in.meta);

end

