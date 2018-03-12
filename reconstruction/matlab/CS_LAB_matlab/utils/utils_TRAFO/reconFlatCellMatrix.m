function reconCell = reconFlatCellMatrix( flatCell, index )
% reconstruct a flattened cell array
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% maximal 4D cell array can be reconstructed
if(nargin == 2)
    if(isempty(index))
        reconCell = flatCell;
        return;
    else
        flatCell = cat(1,flatCell,index);
    end
end

ndim = max(cellfun(@(x) size(x,2), flatCell(2,:)));
dimLevel = max(cell2mat(cellfun(@(x) x(1,:), flatCell(2,:), 'UniformOutput', false).'),[],1);

reconCell = cell(dimLevel);

for i=1:numel(reconCell)
    % get current row, col, slice, page, ...
    [row, col, slice, page] = ind2sub(size(reconCell),i);
    pos = [row, col, slice, page];
    pos = pos(1:ndim);

    currCells = cell(2,size(flatCell,2));
    tmp = false(1,size(flatCell,2));
    currCells_idx = 1;
    for iInner=1:size(flatCell,2)
        if(ismember(flatCell{2,iInner}(1,:),pos,'rows'))
            currCells{1,currCells_idx} = flatCell{1,iInner};
            currCells{2,currCells_idx} = flatCell{2,iInner}(2:end,:);
            tmp(1,currCells_idx) = isempty(currCells{2,currCells_idx});
            currCells_idx = currCells_idx + 1;
        end
    end
    currCells = currCells(:,1:currCells_idx-1);
    tmp = tmp(:,1:currCells_idx-1);
%     currCells_idx = cellfun(@(x) ismember(x(1,:),pos,'rows'), flatCell(2,:));
%     currCells = flatCell(:,currCells_idx);
%     currCells(2,:) = cellfun(@(x) x(2:end,:), currCells(2,:), 'UniformOutput', false);
%     
%     tmp = cellfun(@isempty, currCells(2,:));
    if(any(tmp))
        for j=find(tmp)
            cmd = 'reconCell{pos(1)';
            for k=2:length(pos)
                cmd = sprintf('%s,pos(%d)',cmd,k);
            end
            cmd = sprintf('%s} = currCells{1,j};',cmd);
            eval(cmd); % reconCell{pos} = currCells{1,j};
        end
    else
        helpCell = reconFlatCellMatrix(currCells);
        cmd = 'reconCell{pos(1)';
        for k=2:length(pos)
            cmd = sprintf('%s,pos(%d)',cmd,k);
        end
        cmd = sprintf('%s} = helpCell;',cmd);
        eval(cmd); % reconCell{pos} = helpCell(:);
    end
    
end