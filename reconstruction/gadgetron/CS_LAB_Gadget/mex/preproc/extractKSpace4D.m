function [kSpace, evalMask] = extractKSpace4D(dData, iLC, iEvalInfoMask, dim, inLoops, phaseLCcol)
% extract kSpace from acquired 4D data
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
nTime = dim(4);
iTotalLines = dim(1);
kSpace = cell(1,nCha);
evalMask = cell(1,nCha);

for iCha = 1:nCha
    kSpace{1,iCha} = complex(zeros(iTotalLines, iLC(1, 1), nZ, nTime,'single'), zeros(iTotalLines, iLC(1, 1), nZ, nTime,'single'));
    evalMask{1,iCha} = zeros(iTotalLines, 1, nZ, nTime);
%     fprintf(1, 'Preprocessing Channel %2u\n', iCha);
    dChaData = dData(iCha:nCha:end, :);
    iChaLC   = iLC  (iCha:nCha:end, :);
    if(~isempty(iEvalInfoMask))
        iEvalInfoMaskLC = iEvalInfoMask(iCha:nCha:end, :);
    end
    for t = 1:nTime
        dPhaData = dChaData(iChaLC(:, phaseLCcol) == t - 1, :);
        iPhaLC   = iChaLC  (iChaLC(:, phaseLCcol) == t - 1, :);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMaskPhaLC = iEvalInfoMaskLC(iChaLC(:,phaseLCcol) == t - 1, :);
        end
        
        for iPar = 1:nZ
            dParData = dPhaData(iPhaLC(:, 6) == iPar - 1, :);
            iParLC   = iPhaLC  (iPhaLC(:, 6) == iPar - 1, :);
            if(~isempty(iEvalInfoMask))
                iEvalInfoMaskParLC = iEvalInfoMaskPhaLC(iPhaLC(:,6) == iPar - 1, :);
            end
            for iSeg = 1:nSeg
                dSegData = dParData(iParLC(:,11) == iSeg - 1, :);
                iSegLC = iParLC(iParLC(:,11) == iSeg - 1, :);
                if(~isempty(iEvalInfoMask))
                    iEvalInfoMaskSegLC = iEvalInfoMaskParLC(iParLC(:,11) == iSeg - 1, :);
                end
                for iLine = 1:length(iSegLC(:,3))
                    kSpace{1,iCha}(iSegLC(iLine,3) + 1, :, iPar, t) = dSegData(iLine, :);
                    if(~isempty(iEvalInfoMask) && ~isempty(iEvalInfoMaskSegLC))
                        evalMask{1,iCha}(iSegLC(iLine,3) + 1, 1, iPar, t) = iEvalInfoMaskSegLC(iLine,1);
                    end
                    if(~isempty(iEvalInfoMask) && ~isempty(iEvalInfoMaskParLC) && bitget(iEvalInfoMaskParLC(iLine,1),25)) % switched readout direction
    %                 if(iSeg == 2) % switched readout
                        kSpace{1,iCha}(iSegLC(iLine,3) + 1, :, iPar, t) = kSpace{1,iCha}(iSegLC(iLine,3) + 1, end:-1:1, iPar, t);
    %                     kSpace{1,iCha}(iSegLC(iLine, 3) + 1, :, iPar) = kSpace{1,iCha}(iSegLC(iLine, 3) + 1, end:-1:1, iPar);
                    end
                end
                dispProgress('Extracting K-Space', ((inLoops(1)-1)*inLoops(2)*nCha*nZ*nSeg*nTime + (inLoops(2)-1)*nCha*nZ*nSeg*nTime + (iCha-1)*nZ*nSeg*nTime + (t-1)*nZ*nSeg + iPar*nSeg + iSeg)/inLoops(3));
            end
        end
    end
    dispProgress('Extracting K-Space', ((inLoops(1)-1)*inLoops(2)*nCha*nZ*nSeg*nTime + (inLoops(2)-1)*nCha*nZ*nSeg*nTime + iCha*nZ*nSeg*nTime)/inLoops(3));
end

end
