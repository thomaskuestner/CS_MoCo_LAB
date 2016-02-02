function out = ones( in, sPrecision )
%ONES for initialization of nested cell arrays
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

out = cellfun(@(x) ones(x, sPrecision), in.data{2,1}, 'UniformOutput', false);
out = reconFlatCellMatrix(out, in.data{1,1});
out = TRAFO(out,in.meta);

end

