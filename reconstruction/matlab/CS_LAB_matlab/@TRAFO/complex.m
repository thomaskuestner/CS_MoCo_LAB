function out = complex( inA, inB )
% complex for initialization of nested cell arrays
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

[outA, idx] = flattenCellMatrix(inA.data);
outB = flattenCellMatrix(inB.data);
out = cellfun(@(x,y) complex(x,y), outA, outB, 'UniformOutput', false);
out = reconFlatCellMatrix(out,idx);
out = TRAFO(out,inA.meta);

end

