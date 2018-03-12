function out = size( in )
%SIZE for nested cell arrays
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

if(nargin > 1)
    error('TRAFO::size: size of specific dimension cannot be extracted');
end

[flatCell, idx] = flattenCellMatrix(in.data);
out = cell(2,1);
out{1,1} = idx;
% out{2,1} = cell(1,size(flatCell,2));
% for i=1:size(flatCell,2)
%     out{2,1}(1,i) = size(flatCell{1,i});
% end
out{2,1} = cellfun(@(x) size(x), flatCell, 'UniformOutput', false);
out = TRAFO(out,in.meta);

end

