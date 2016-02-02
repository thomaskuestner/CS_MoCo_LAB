function [out,pos] = max( in, plonk, dim )
%MAX for nested cells
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

% return scalar value
[out, ~] = flattenCellMatrix(in.data);
if(nargin == 1)
    out = max(cellfun(@(x) max(x(:)), out));
else
    [out, pos] = cellfun(@(x) max(x,[],dim), out, 'UniformOutput', false);
end

end

