function out = sum( in, dim )
%SUM for nested cells
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

if(nargin < 2)
    error('TRAFO::sum: Use plus instead');
end

if(dim == 1 || dim == 2 || dim == 3) % summation inside cell
    [out, idx] = flattenCellMatrix(in.data);
    out = cellfun(@(x) sum(x,dim), out, 'UniformOutput',false);
    out = reconFlatCellMatrix(out,idx);
    out = TRAFO(out, in.meta);
    
if(dim == 4)
    % summation along channels
    data = in.data;
    out = zeros(size(in));
    [out, idx] = flattenCellMatrix(out.data(1));
    for i=2:length(data)
        tmp = flattenCellMatrix(in.data(i));
        out = cellfun(@(x,y) x + y, out, tmp, 'UniformOutput', false);
    end
    out = reconFlatCellMatrix(out,idx);
%     if iscell(in.meta)
%         out = TRAFO(out, in.meta(1));
%     else 
        out = TRAFO(out, in.meta);
%     end
end

end

