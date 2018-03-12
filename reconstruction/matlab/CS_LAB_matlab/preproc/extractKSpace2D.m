function [kSpace, evalMask] = extractKSpace2D(dData, iLC, iEvalInfoMask, dim, inLoops)
% extract kSpace from acquired 2D data
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% Check the input data
% These loop counters have to be constant in one slice
% iSingularLCs = iLC(:,[1, 2, 5, 6, 7, 9, 10, 11]); % Sam, Cha, Sli, Par, Eco, Rep, Set, Seg
% if any(min(iSingularLCs) ~= max(iSingularLCs))
%     error('CS_reconstruction()::extractKSpace2D()','Input data contains data from multiple sample sizes/channel counts/slices/partitions/echos/repetitions/sets/segments!');
% end

% extract kSpace
nCha = dim(5);
nTime = dim(4);
nSeg = double(max(iLC(:,11)) - min(iLC(:,11)) + 1);
iTotalLines = dim(1);
kSpace = cell(1,nCha);
evalMask = cell(1,nCha);

for iCha = 1:nCha
    kSpace{iCha} = complex(zeros(iTotalLines, iLC(1, 1), nTime,'single'), zeros(iTotalLines, iLC(1, 1), nTime,'single'));
    evalMask{iCha} = zeros(iTotalLines, 1, nTime);
    
%     fprintf(1, 'Preprocessing Channel %2u\n', iCha);
    dChaData = dData(iCha:nCha:end, :);
    iChaLC   = iLC  (iCha:nCha:end, :);
    if(~isempty(iEvalInfoMask))
        iEvalInfoMaskLC = iEvalInfoMask(iCha:nCha:end, :);
    end
    for t = 1:nTime
        dPhaData = dChaData(iChaLC(:, 8) == t - 1, :);
        iPhaLC   = iChaLC  (iChaLC(:, 8) == t - 1, :);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMaskPhaLC = iEvalInfoMaskLC(iChaLC(:,8) == t - 1, :);
        end
        for iSeg = 1:nSeg
            dSegData = dPhaData(iPhaLC(:,11) == iSeg - 1, :);
            iSegLC = iPhaLC(iPhaLC(:,11) == iSeg - 1, :);
            if(~isempty(iEvalInfoMask))
                iEvalInfoMaskSegLC = iEvalInfoMaskPhaLC(iPhaLC(:,11) == iSeg - 1, :);
            end        
            for iLine = 1:length(iSegLC(:,3))
%                 if(~isempty(iEvalInfoMask) && bitget(iEvalInfoMaskSegLC(iLine,1),22))
%                     continue; % reference lines in EPI
% % %                     iSegIdx = find(iSegLC(:,3) == iLine - 1);
% %                     tmpIdx = iSegLC(iLine,3) + 1;
% %                     lIdx = false(length(iSegLC(iLine,3)),1);
% %                     for i=1:length(iSegLC(iLine,3))
% %                         lIdx(i) = ~bitget(iEvalInfoMaskSegLC(iSegIdx(i),1),22); % reference lines in EPI
% %                     end
% %                     kIdx = tmpIdx(lIdx);
% %                     iSegLine = 1;
%                 end
                kSpace{iCha}(iSegLC(iLine,3) + 1, :, t) = dSegData(iLine, :); 
                if(~isempty(iEvalInfoMask) && ~isempty(iEvalInfoMaskSegLC))
                    evalMask{iCha}(iSegLC(iLine,3) + 1, 1, t) = iEvalInfoMaskSegLC(iLine,1);
                end
                if(~isempty(iEvalInfoMask) && bitget(iEvalInfoMaskSegLC(iLine,1),25)) % switched readout direction
%                 if(iSeg == 2) % switched readout
                    kSpace{iCha}(iSegLC(iLine,3) + 1, :, t) = kSpace{iCha}(iSegLC(iLine,3) + 1, end:-1:1, t);
                end
            end
            dispProgress('Extracting K-Space', ((inLoops(1)-1)*inLoops(5)*inLoops(6)*nCha*nTime*nSeg + (inLoops(2)-1)*inLoops(6)*nCha*nTime*nSeg + (inLoops(3)-1)*nCha*nTime*nSeg + (iCha-1)*nTime*nSeg + t*nSeg + iSeg)/inLoops(4));
        end
%         dKSpace(:,:,t) = ifftshift(dKSpace(:,:,t)); % change quadrants
    end
%     dKSpace(1:2:end, :, :) = - dKSpace(1:2:end, :, :);
%     dKSpace(:, 1:2:end, :) = - dKSpace(:, 1:2:end, :); % generates a checker board with sign changes
    
    dispProgress('Extracting K-Space', ((inLoops(1)-1)*inLoops(5)*inLoops(6)*nCha*nTime*nSeg + (inLoops(2)-1)*inLoops(6)*nCha*nTime*nSeg + (inLoops(3)-1)*nCha*nTime*nSeg + iCha*nTime*nSeg)/inLoops(4));
end

end

