function [kSpace, evalMask] = extractKSpace3D(dData, iLC, iEvalInfoMask, dim, inLoops)
% extract kSpace from acquired 3D data
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% Check the input data
% These loop counters have to be constant in one slice
% iSingularLCs = iLC(:,[1, 2, 5, 7, 9, 10, 11]); % Sam, Cha, Sli, Eco, Rep, Set, Seg
% if any(min(iSingularLCs) ~= max(iSingularLCs))
%     error('CS_reconstruction()::extractKSpace3D()','Input data contains data from multiple sample sizes/channel counts/slices/echos/repetitions/sets/segments!');
% end

% extract kSpace
nCha = dim(5);
nZ = dim(3);
nSeg = double(max(iLC(:,11)) - min(iLC(:,11)) + 1);
iTotalLines = dim(1);
kSpace = cell(1,nCha);
evalMask = cell(1,nCha);

for iCha = 1:nCha
    kSpace{1,iCha} = (complex(zeros(iTotalLines, iLC(1, 1), nZ,'single'), zeros(iTotalLines, iLC(1, 1), nZ,'single')));
    evalMask{1,iCha} = zeros(iTotalLines, 1, nZ);
%     fprintf(1, 'Preprocessing Channel %2u\n', iCha);
    dChaData = dData(iCha:nCha:end, :);
    iChaLC   = iLC  (iCha:nCha:end, :);
    if(~isempty(iEvalInfoMask))
        iEvalInfoMaskLC = iEvalInfoMask(iCha:nCha:end, :);
    end
    for iPar = 1:nZ
        dParData = dChaData(iChaLC(:, 6) == iPar - 1, :);
        iParLC   = iChaLC  (iChaLC(:, 6) == iPar - 1, :);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMaskParLC = iEvalInfoMaskLC(iChaLC(:,6) == iPar - 1, :);
        end
        for iSeg = 1:nSeg
            dSegData = dParData(iParLC(:,11) == iSeg - 1, :);
            iSegLC = iParLC(iParLC(:,11) == iSeg - 1, :);
            if(~isempty(iEvalInfoMask))
                iEvalInfoMaskSegLC = iEvalInfoMaskParLC(iParLC(:,11) == iSeg - 1, :);
            end
            for iLine = 1:length(iSegLC(:,3))
                kSpace{1,iCha}(iSegLC(iLine,3) + 1, :, iPar) = dSegData(iLine, :);
                if(~isempty(iEvalInfoMask) && ~isempty(iEvalInfoMaskSegLC))
                    evalMask{1,iCha}(iSegLC(iLine,3) + 1, 1, iPar) = iEvalInfoMaskSegLC(iLine,1);
                end
                if(~isempty(iEvalInfoMask) && ~isempty(iEvalInfoMaskParLC) && bitget(iEvalInfoMaskParLC(iLine,1),25)) % switched readout direction
%                 if(iSeg == 2) % switched readout
                    kSpace{1,iCha}(iSegLC(iLine,3) + 1, :, iPar) = kSpace{1,iCha}(iSegLC(iLine,3) + 1, end:-1:1, iPar);
%                     kSpace{1,iCha}(iSegLC(iLine, 3) + 1, :, iPar) = kSpace{1,iCha}(iSegLC(iLine, 3) + 1, end:-1:1, iPar);
                end
            end
            dispProgress('Extracting K-Space', ((inLoops(1)-1)*inLoops(2)*nCha*nZ*nSeg + (inLoops(2)-1)*nCha*nZ*nSeg + (iCha-1)*nZ*nSeg + iPar*nSeg + iSeg)/inLoops(3));
        end
    end
    dispProgress('Extracting K-Space', ((inLoops(1)-1)*inLoops(2)*nCha*nZ*nSeg + (inLoops(2)-1)*nCha*nZ*nSeg + iCha*nZ*nSeg)/inLoops(3));
end

end
