function [kSpace, evalMask] = extractKSpace_main(dData, iLC, iEvalInfoMask, measPara)
% extract kSpace data
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% set repetitions
if(measPara.LCall(4) ~= 1)
    iRepetitions = 1:measPara.LCall(4);
else
    iRepetitions = 1;
end

% set averages
if(measPara.LCall(2) ~= 1)
    iAverages = 1:measPara.LCall(2);
else
    iAverages = 1;
end

nSeg = double(max(iLC(:,11)) - min(iLC(:,11)) + 1);

fprintf(1, 'Preprocessing Channels\n');    
if(strcmp(measPara.dimension,'2D'))
    % cell: slice - cha
    % in cell: k_y - k_x - (t)
    if(measPara.LCall(1) == 1)
        % one 2D slice was acquired
        iSlices = 1;
    else
        % several 2D slices were acquired        
        iSlices = 1:measPara.LCall(1);
    end

    kSpace = cell(length(iSlices),measPara.dim(5),length(iRepetitions),length(iAverages));
    evalMask = cell(length(iSlices),measPara.dim(5),length(iRepetitions),length(iAverages));
    
    % extract phase correction data
    if(~isempty(iEvalInfoMask))
        lPhaseCorr = false(size(iLC,1),1);
        for iPhaseCorr=1:length(lPhaseCorr)
            lPhaseCorr(iPhaseCorr) = bitget(iEvalInfoMask(iPhaseCorr,1),22);
        end
        iLC = iLC(~lPhaseCorr,:);
        dData = dData(~lPhaseCorr,:);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMask = iEvalInfoMask(~lPhaseCorr,:);
        end
    end

    totalCount = length(iRepetitions)*length(iAverages)*measPara.dim(5)*length(iSlices)*measPara.dim(4)*nSeg;
    dispProgress('Extracting K-Space', 0, totalCount);
    for iRep=1:length(iRepetitions)
        dDataRep = dData(iLC(:,9) == iRepetitions(iRep)-1, :);
        iLCRep = iLC(iLC(:,9) == iRepetitions(iRep)-1, :);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMaskRepLC = iEvalInfoMask(iLC(:,9) == iRepetitions(iRep)-1, :);
        end
        for iAvg=1:length(iAverages)
            dDataAvg = dDataRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
            iLCAvg = iLCRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
            if(~isempty(iEvalInfoMask))
                iEvalInfoMaskAvgLC = iEvalInfoMaskRepLC(iLCRep(:,4) == iAverages(iAvg)-1, :);
            end
            for iSli=1:length(iSlices)
                dDataSli = dDataAvg(iLCAvg(:,5) == iSlices(iSli)-1, :);
                iLCSli = iLCAvg(iLCAvg(:,5) == iSlices(iSli)-1, :);
                if(~isempty(iEvalInfoMask))
                    iEvalInfoMaskSliLC = iEvalInfoMaskAvgLC(iLCAvg(:,5) == iSlices(iSli)-1, :);
                else
                    iEvalInfoMaskSliLC = [];
                end

                [kSpace(iSli,:,iRep,iAvg), evalMask(iSli,:,iRep,iAvg)] = extractKSpace2D(dDataSli, iLCSli, iEvalInfoMaskSliLC, measPara.dim, [iRep, iAvg, iSli, totalCount, length(iAverages), length(iSlices)]);
                dispProgress('Extracting K-Space', ((iRep-1)*length(iAverages)*measPara.dim(5)*length(iSlices)*measPara.dim(4)*nSeg + (iAvg-1)*measPara.dim(5)*length(iSlices)*measPara.dim(4)*nSeg + iSli*measPara.dim(5)*measPara.dim(4)*nSeg)/totalCount);
            end
            dispProgress('Extracting K-Space', ((iRep-1)*length(iAverages)*measPara.dim(5)*length(iSlices)*measPara.dim(4)*nSeg + iAvg*measPara.dim(5)*length(iSlices)*measPara.dim(4)*nSeg)/totalCount);
        end
        dispProgress('Extracting K-Space', iRep/length(iRepetitions));
    end
    dispProgress('Extracting K-Space', 'Close');

elseif(strcmp(measPara.dimension,'3D'))
    % cell: 1 - cha
    % in cell: k_y - k_x - k_z
    
    kSpace = cell(1,measPara.dim(5),length(iRepetitions),length(iAverages));
    evalMask = cell(1,measPara.dim(5),length(iRepetitions),length(iAverages));
    totalCount = length(iRepetitions)*measPara.dim(5)*measPara.dim(3)*nSeg;
    dispProgress('Extracting K-Space', 0, totalCount);
    for iRep=1:length(iRepetitions)
        dDataRep = dData(iLC(:,9) == iRepetitions(iRep)-1, :);
        iLCRep = iLC(iLC(:,9) == iRepetitions(iRep)-1, :);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMaskRepLC = iEvalInfoMask(iLC(:,9) == iRepetitions(iRep)-1, :);
        end
        for iAvg=1:length(iAverages)
            dDataAvg = dDataRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
            iLCAvg = iLCRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
            if(~isempty(iEvalInfoMask))
                iEvalInfoMaskAvgLC = iEvalInfoMaskRepLC(iLCRep(:,4) == iAverages(iAvg)-1, :);
            else
                iEvalInfoMaskAvgLC = [];
            end
            [kSpace(1,:,iRep,iAvg), evalMask(1,:,iRep,iAvg)] = extractKSpace3D(dDataAvg, iLCAvg, iEvalInfoMaskAvgLC, measPara.dim, [iRep,iAvg,totalCount]);
            dispProgress('Extracting K-Space', ((iRep-1)*length(iAverages) + iAvg)/(length(iRepetitions)*length(iAverages)));
        end
        dispProgress('Extracting K-Space', iRep/length(iRepetitions));
    end
    dispProgress('Extracting K-Space', 'Close');

elseif(strcmp(measPara.dimension,'4D'))
    % cell: 1 - cha
    % in cell: k_y - k_x - k_z - t
    
    % extract prescans/Navi for CS_Trufi
    if(max(iLC(:,1)) ~= min(iLC(:,1)))
        dData = dData(iLC(:,1) == measPara.dim(2),:);
        iLC = iLC(iLC(:,1) == measPara.dim(2),:);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMask = iEvalInfoMask(iLC(:,1) == measPara.dim(2),:);
        end
    end
    if(strcmp(measPara.sequenceName,'CS_Trufi'))
        phaseLCcol = 10;
    else
        phaseLCcol = 8;
    end
    kSpace = cell(1,measPara.dim(5),length(iRepetitions),length(iAverages));
    evalMask = cell(1,measPara.dim(5),length(iRepetitions),length(iAverages));
    totalCount = length(iRepetitions)*measPara.dim(5)*measPara.dim(3)*nSeg*measPara.dim(4);
    dispProgress('Extracting K-Space', 0, totalCount);
    for iRep=1:length(iRepetitions)
        dDataRep = dData(iLC(:,9) == iRepetitions(iRep)-1, :);
        iLCRep = iLC(iLC(:,9) == iRepetitions(iRep)-1, :);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMaskRepLC = iEvalInfoMask(iLC(:,9) == iRepetitions(iRep)-1, :);
        end
        for iAvg=1:length(iAverages)
            dDataAvg = dDataRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
            iLCAvg = iLCRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
            if(~isempty(iEvalInfoMask))
                iEvalInfoMaskAvgLC = iEvalInfoMaskRepLC(iLCRep(:,4) == iAverages(iAvg)-1, :);
            else
                iEvalInfoMaskAvgLC = [];
            end
            [kSpace(1,:,iRep,iAvg), evalMask(1,:,iRep,iAvg)] = extractKSpace4D(dDataAvg, iLCAvg, iEvalInfoMaskAvgLC, measPara.dim, [iRep,iAvg,totalCount], phaseLCcol);
            dispProgress('Extracting K-Space', ((iRep-1)*length(iAverages) + iAvg)/(length(iRepetitions)*length(iAverages)));
        end
        dispProgress('Extracting K-Space', iRep/length(iRepetitions));
    end
    dispProgress('Extracting K-Space', 'Close');

end

end

