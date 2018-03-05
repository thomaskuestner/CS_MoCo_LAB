function [flatArray, index] = flattenCellMatrix(nestedArray)
% flatten nested (TRAFO) cell array
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------
    
    if(nnz(size(nestedArray) ~= 1) > 1 || any(any(any(any(cellfun(@iscell,nestedArray))))))
        [flatArray,maxdim] = flattenCell(nestedArray);
    
        index = flatArray(2,:);
        index = cellfun(@(x) x(:,1:maxdim), index, 'UniformOutput', false);
        flatArray = flatArray(1,:);
    else
        index = [];
        flatArray = nestedArray(:).';
    end
    
end

function [flatArray,maxdim] = flattenCell(nestedArray,cnt)

% maximal 4D cell array as input

narginchk(1,2);
if(~exist('cnt','var'))
    cnt = 1;
end
if ~iscell(nestedArray),
    error('Must be a cell array.');
end
flatArray{2,1} = [];
for i=1:numel(nestedArray)
    if iscell(nestedArray{i})
        [y, maxdim] = flattenCell(nestedArray{i},cnt+1);
        [flatArray{:,end+1:end+size(y,2)}] = deal(y{:});
        for j=size(flatArray,2):-1:size(flatArray,2)-size(y,2)+1
            [row, col, slice, page] = ind2sub(size(nestedArray),i);
            flatArray{2,j}(cnt,:) = [row, col, slice, page];
            maxdim = max([maxdim, ndims(nestedArray)]);
        end
    else
        flatArray{1,end+1} = nestedArray{i};
        flatArray{2,end} = zeros(cnt,4);
        maxdim = ndims(nestedArray);
        [row, col, slice, page] = ind2sub(size(nestedArray),i);
        flatArray{2,end}(cnt,:) = [row, col, slice, page];
    end
end
flatArray(:,1) = [];

end

