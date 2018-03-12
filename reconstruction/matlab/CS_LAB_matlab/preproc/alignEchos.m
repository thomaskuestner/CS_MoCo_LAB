function kSpace = alignEchos(kSpace, evalMask, measPara, dData, iLC, iEvalInfoMask)
% align echos in EPI sequences
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% [dData, iLC, iEvalInfoMask] = fMeasRead(fullfile(path,[filename,'.dat']), 'Set', 0, 'Seg', 1, 'Smp', [21 inf]); % old code

if(isempty(dData) || isempty(iLC))
    return;
end

% extract phase correction data
if(~isempty(iEvalInfoMask))
    lPhaseCorr = false(size(iLC,1),1);
    for iPhaseCorr=1:length(lPhaseCorr)
        lPhaseCorr(iPhaseCorr) = bitget(iEvalInfoMask(iPhaseCorr,1),22);
    end
    iLC = iLC(lPhaseCorr,:);
    dData = dData(lPhaseCorr,:);
    if(~isempty(iEvalInfoMask))
        iEvalInfoMask = iEvalInfoMask(lPhaseCorr,:);
    end
end

% set repetitions
if(measPara.LCall(4) ~= 1)
    iRepetitions = 1:measPara.LCall(4);
else
    iRepetitions = 1;
end

% set averages for reference scans
iAverages = 1:double(max(iLC(:,4)) - min(iLC(:,4)) + 1); % for reference scans just always two averages (0 0 1)

% set averages
if(measPara.LCall(2) ~= 1)
    iAveragesData = 1:measPara.LCall(2);
else
    iAveragesData = 1;
end

nSeg = double(max(iLC(:,11)) - min(iLC(:,11)) + 1);
iSegRef = [1;0;1];
iAvgRef = [0;0;1];

fprintf(1, 'Aligning Echos\n');    
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

%     phase = cell(length(iSlices),measPara.dim(5),length(iRepetitions));
%     xVec = -ceil(measPara.dim(2)/2):ceil(measPara.dim(2)/2)-1;

    

    nCha = measPara.dim(5);
    totalCount = length(iRepetitions)*length(iAveragesData)*measPara.dim(5)*length(iSlices)*measPara.dim(4);
    dispProgress('Phase Correction', 0, totalCount);
    for iRep=1:length(iRepetitions)
        dDataRep = dData(iLC(:,9) == iRepetitions(iRep)-1, :);
        iLCRep = iLC(iLC(:,9) == iRepetitions(iRep)-1, :);
        if(~isempty(iEvalInfoMask))
            iEvalInfoMaskRepLC = iEvalInfoMask(iLC(:,9) == iRepetitions(iRep)-1, :);
        end
%         for iAvg=1:length(iAverages)
%             dDataAvg = dDataRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
%             iLCAvg = iLCRep(iLCRep(:,4) == iAverages(iAvg)-1, :);
            
        for iAvgData = 1:length(iAveragesData)
            
            for iSli=1:length(iSlices)
                dDataSli = dDataRep(iLCRep(:,5) == iSlices(iSli)-1, :);
                iLCSli = iLCRep(iLCRep(:,5) == iSlices(iSli)-1, :);
                if(~isempty(iEvalInfoMask))
                    iEvalInfoMaskSliLC = iEvalInfoMaskRepLC(iLCRep(:,5) == iSlices(iSli)-1, :);
                end
                for iCha = 1:nCha
                    phaseCorr = complex(zeros(3, iLCSli(1, 1), measPara.dim(4)), zeros(3, iLCSli(1, 1), measPara.dim(4)));
                    dChaData = dDataSli(iCha:nCha:end, :);
                    iChaLC   = iLCSli(iCha:nCha:end, :);
                    if(~isempty(iEvalInfoMask))
                        iEvalInfoMaskLC = iEvalInfoMaskSliLC(iCha:nCha:end, :);
                    end
                    for t = 1:measPara.dim(4)
                        readoutDir = zeros(3,1);
                        dPhaData = dChaData(iChaLC(:, 8) == t - 1, :);
                        iPhaLC   = iChaLC  (iChaLC(:, 8) == t - 1, :);
                        if(~isempty(iEvalInfoMask))
                            iEvalInfoMaskPhaLC = iEvalInfoMaskLC(iChaLC(:,8) == t - 1, :);
                        end
                        % find inverse readouts
                        lMask = abs(kSpace{iSli,iCha,iRep,iAvgData}(:,:,t)) > 0;
                        lMaskInverse = false(size(lMask));
                        for iLine = 1:measPara.dim(1)
                            if(~isempty(evalMask{iSli,iCha,iRep,iAvgData}) && bitget(evalMask{iSli,iCha,iRep,iAvgData}(iLine,1,t),25)) % switched readout direction
%                                     if(~isempty(evalMask{iSli,iCha,iRep,iAvgData}) && iSeg == 1)
                                lMaskInverse(iLine,:) = lMask(iLine,:) & true;
                            end
                        end
                        
                        % cut out reference scan corresponding to average
                        % loop of data
                        dAvgDataData = dPhaData((iAvgData-1)*length(iAvgRef)+1:iAvgData*length(iAvgRef), :);
                        iAvgDataLC = iPhaLC((iAvgData-1)*length(iAvgRef)+1:iAvgData*length(iAvgRef), :);
                        if(~isempty(iEvalInfoMask))
                            iEvalInfoMaskAvgDataLC = iEvalInfoMaskPhaLC((iAvgData-1)*length(iAvgRef)+1:iAvgData*length(iAvgRef), :);
                        end
                        
                        for iAvg = 1:length(iAverages)
                            dAvgData = dAvgDataData(iAvgDataLC(:,4) == iAverages(iAvg)-1, :);
                            iAvgLC   = iAvgDataLC(iAvgDataLC(:,4) == iAverages(iAvg)-1, :);
                            if(~isempty(iEvalInfoMask))
                                iEvalInfoMaskAvgLC = iEvalInfoMaskAvgDataLC(iAvgDataLC(:,4) == iAverages(iAvg)-1, :);
                            end                            
                            for iSeg = 1:nSeg
                                dSegData = dAvgData(iAvgLC(:,11) == iSeg - 1, :);
                                iSegLC = iAvgLC(iAvgLC(:,11) == iSeg - 1, :);
                                if(~isempty(iEvalInfoMask))
                                    iEvalInfoMaskSegLC = iEvalInfoMaskAvgLC(iAvgLC(:,11) == iSeg - 1, :);
                                end
                                for iLine = 1:length(iSegLC(:,3))
            %                         phaseCorr(iPhaLC(iLine, 3) + 1, :, t) = dPhaData(iLine, :); % eventuell anpassen für drei echo lines mit selbem iPhaLC
            %                         readoutDir(iPhaLC(iLine,3) + 1, 1) = 1;
                                    if(~isempty(iEvalInfoMask) && bitget(iEvalInfoMaskSegLC(iLine,1),22))
%                                         tmpIdx = iSegLC(iLine,3) + 1;
%                                         lIdx = false(length(iSegLC(iLine,3)),1);
%                                         for i=1:length(iSegLC(iLine,3))
%                                             lIdx(i) = bitget(iEvalInfoMaskSegLC(iSegIdx(i),1),22); % reference lines in EPI
%                                         end
%                                         kIdx = tmpIdx(lIdx);
%                                         refData = dSegData(lIdx, :);
%                                         refData = refData(iAvgData, :);
                                        phaseCorr(iSegRef == iSeg - 1 & iAvgRef == iAvg - 1, :, t) = dSegData(iLine, :);
                                        readoutDir(iSegRef == iSeg - 1 & iAvgRef == iAvg - 1, 1) = 1;
                                        if(~isempty(iEvalInfoMask) && bitget(iEvalInfoMaskSegLC(iLine,1),25)) % switched readout direction
%                                         if(iSeg == 2) % switched readout direction
                                            phaseCorr(iSegRef == iSeg - 1 & iAvgRef == iAvg - 1, :, t) = phaseCorr(iSegRef == iSeg - 1 & iAvgRef == iAvg - 1, end:-1:1, t);
                                            readoutDir(iSegRef == iSeg - 1 & iAvgRef == iAvg - 1, 1) = -1;
                                        end
                                    end

%                                     phaseCorr(iLine, :, t) = dSegData(iLine, :);
%                                     readoutDir(iLine, 1) = 1;
    %                                 if(~isempty(iEvalInfoMask) && bitget(iEvalInfoMaskPhaLC(iLine,1),25)) % switched readout direction
            %                             phaseCorr(iPhaLC(iLine, 3) + 1, :, t) = phaseCorr(iPhaLC(iLine, 3) + 1, end:-1:1, t);
%                                         readoutDir(iLine, 1) = -1;
%                                     end
                                end

                                
                            end
                        end

                        phasePos = repmat(angle(ifftnshift(mean(phaseCorr(readoutDir == 1,:,t),1),2)),[measPara.dim(1), 1]);
                        phaseNeg = repmat(angle(ifftnshift(mean(phaseCorr(readoutDir == -1,:,t),1),2)),[measPara.dim(1), 1]); % mean just to ensure to have one single line

                        hybridImg = ifftnshift(kSpace{iSli,iCha,iRep,iAvgData}(:,:,t),2);
                        kSpace{iSli,iCha,iRep,iAvgData}(:,:,t) = fftnshift(abs(hybridImg) .* exp(1i * (angle(hybridImg) - phasePos .* ~lMaskInverse - phaseNeg .* lMaskInverse)),2);
    %                     kSpace{iSli,iCha,iRep}(:,:,t) = fftnshift(abs(hybridImg) .* exp(1i * (angle(hybridImg) - phasePos .* ~lMaskInverse - phaseNeg .* lMaskInverse)),2);

    %                     phaseDelta = zeros(1,measPara.dim(2));
    %                     rowIdx = find(sum(phaseCorr(:,:,t),2),1,'first');

    % %                     % method a): maximum correction
    % %                     [~, maxRef] = max(phaseCorr(rowIdx,:,t));
    % %                     [~, maxImg] = max(kSpace{iSli,iCha,iRep}(rowIdx,:,t));
    % %                    
    % %                     diffPos = maxImg - maxRef;
    % %                     kSpace{iSli,iCha,iRep}(:,:,t) = circshift(kSpace{iSli,iCha,iRep}(:,:,t), [0, diffPos]) .* lMaskInverse + kSpace{iSli,iCha,iRep}(:,:,t) .* ~lMaskInverse;

    %                     echoRefL = fftnshift(zpad(echoRef,[1, 2*size(echoRef,2)]),2);
    %                     echoImgL = fftnshift(zpad(echoImg,[1, 2*size(echoImg,2)]),2);

                        % Conclusion:
                        % - Method b) better than a)
                        % - correct for inverse readout! (end:-1:1) needed
                        % - alignEchos needed!

                        % method b): linear phase correction
    %                     echoRef = ifftnshift(phaseCorr(rowIdx,:,t),2);
    %                     echoImg = ifftnshift(kSpace{iSli,iCha,iRep}(rowIdx,:,t),2);  
    %                     % just take image and neglect noise part
    %                     blocks = find_blocks(abs(echoRef) > median(abs(echoRef)) & abs(echoImg) > median(abs(echoImg))); % & [abs(diff(angle(echoRef))) < median(diff(angle(echoRef))),false] & [abs(diff(angle(echoImg))) < median(diff(angle(echoImg))), false]);
    %                     blocks = blocks(:,abs(blocks(1,:) - blocks(2,:)) > 1 & blocks(1,:) <= ceil(measPara.dim(2)/2) & blocks(2,:) >= ceil(measPara.dim(2)/2));
    %                     colIdx = blocks(1):blocks(2); % rough estimate for inner part
    %                     
    %                     % estimate phase difference
    %                     phaseRef = angle(echoRef(1,colIdx));
    %                     phaseImg = angle(echoImg(1,colIdx));
    %                     
    %                     % estimate linear phase
    %                     A = [ones(length(xVec(colIdx)),1) xVec(colIdx).'];
    %                     theta(:,1) = A\(phaseRef.');
    %                     theta(:,2) = A\(phaseImg.');
    %                     phaseDelta(1,colIdx) = A * (theta(:,1) - theta(:,2));
    %                     hybridImg = ifftnshift(kSpace{iSli,iCha,iRep}(:,:,t),2);
    %                     
    %                     phaseDelta = repmat(phaseDelta,[size(hybridImg,1),1]);
    %                     kSpace{iSli,iCha,iRep}(:,:,t) = fftnshift(abs(hybridImg) .* exp(1i * (angle(hybridImg) + phaseDelta .* lMaskInverse)),2);


    % %                     
    % %                     % unwrap phases
    % %                     a = unwrap(phaseRef);
    % %                     
    % %                     
    % %                     figure;
    % %                     subplot(2,2,1); plot(angle(echoRef)); title('phaseRef');
    % %                     subplot(2,2,3); plot(angle(echoImg)); title('phaseImg');
    % %                     subplot(2,2,2); plot(unwrap(angle(echoRef))); title('phaseRef_unwrap');
    % %                     subplot(2,2,4); plot(unwrap(angle(echoImg))); title('phaseImg_unwrap');

                        % ifft -> zeropad -> fft -> find max -> delay -> ifft
                        % -> crop -> fft

    %                     for i=1:length(xVec)
    %     %                     theta = [1 xVec(i); 1 xVec(i)]\[phaseCorr(i); phaseRef(i)];
    %                         theta = pinv([1 xVec(i); 1 xVec(i)]) * [phaseRef(i); phaseImg(i)];
    %                         phaseDelta(i) = theta.' * [1; xVec(i)];
    %                     end

                        dispProgress('Phase Correction', ((iRep-1)*length(iAveragesData)*length(iSlices)*nCha*measPara.dim(4) + (iAvgData-1)*length(iSlices)*nCha*measPara.dim(4) + (iSli-1)*nCha*measPara.dim(4) + (iCha-1)*measPara.dim(4) + t)/totalCount);
                    end
                    dispProgress('Phase Correction', ((iRep-1)*length(iAveragesData)*length(iSlices)*nCha*measPara.dim(4) + (iAvgData-1)*length(iSlices)*nCha*measPara.dim(4) + (iSli-1)*nCha*measPara.dim(4) + iCha*measPara.dim(4))/totalCount);
                end
                dispProgress('Phase Correction', ((iRep-1)*length(iAveragesData)*measPara.dim(5)*length(iSlices)*measPara.dim(4) + (iAvgData-1)*length(iSlices)*nCha*measPara.dim(4) + iSli*measPara.dim(5)*measPara.dim(4))/totalCount);
            end
            dispProgress('Phase Correction', ((iRep-1)*length(iAveragesData)*measPara.dim(5)*length(iSlices)*measPara.dim(4) + iAvgData*length(iSlices)*nCha*measPara.dim(4))/totalCount);
        end
        dispProgress('Phase Correction', iRep/length(iRepetitions));
    end
    dispProgress('Phase Correction', 'Close');
        
elseif(strcmp(measPara.dimension,'3D'))
    error('Not yet supported');
end

end