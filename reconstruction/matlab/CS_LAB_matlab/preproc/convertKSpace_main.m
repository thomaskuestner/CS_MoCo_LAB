function [ kSpace, measPara ] = convertKSpace_main( data )
% convert from ISMRMD to processable k-space format
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

iFreq = unique(double(data.head.number_of_samples));
iY = unique(double(data.head.idx.kspace_encode_step_1));
iZ = unique(double(data.head.idx.kspace_encode_step_2));
iPha = unique(double(data.head.idx.phase));
iSli = unique(double(data.head.idx.slice));
iRep = unique(double(data.head.idx.repetition));
iAvg = unique(double(data.head.idx.average));
iEcho = unique(double(data.head.idx.contrast));
iCha = unique(double(data.head.active_channels));
measPara.dim = [length(iY), iFreq, length(iZ), length(iPha), iCha];
measPara.LCall = [length(iSli), length(iAvg), length(iEcho), length(iRep)];

if(measPara.dim(3) == 1)
    measPara.dimension = '2D';
else
    if(measPara.dim(4) == 1)
        measPara.dimension = '3D';
    else
        measPara.dimension = '4D';
    end    
end

% begin conversion        
kSpace = cell(measPara.LCall(1), measPara.dim(5), measPara.LCall(4), measPara.LCall(2));
if(strcmp(measPara.dimension,'2D'))
    [kSpace{:}] = deal(complex(zeros(measPara.dim(1),iFreq,measPara.dim(4),'single'),zeros(measPara.dim(1),iFreq,measPara.dim(4),'single')));
elseif(strcmp(measPara.dimension,'3D'))
    [kSpace{:}] = deal(complex(zeros(measPara.dim(1),iFreq,measPara.dim(3),'single'),zeros(measPara.dim(1),iFreq,measPara.dim(3),'single')));
elseif(strcmp(measPara.dimension,'4D'))
    [kSpace{:}] = deal(complex(zeros(measPara.dim(1),iFreq,measPara.dim(3),measPara.dim(4),'single'),zeros(measPara.dim(1),iFreq,measPara.dim(3),measPara.dim(4),'single')));
end

for idxRep = 1:measPara.LCall(4)
    for idxAvg = 1:measPara.LCall(2)
        
        if(strcmp(measPara.dimension,'2D')) 
            for idxSli = 1:measPara.LCall(1)
                for idxPha = 1:measPara.dim(4)
                    for idxY = 1:measPara.dim(1)
                        lIdx = double(data.head.idx.repetition) == iRep(idxRep) & ...
                            double(data.head.idx.average) == iAvg(idxAvg) & ...
                            double(data.head.idx.slice) == iSli(idxSli) & ...
                            double(data.head.idx.phase) == iPha(idxPha) & ...
                            double(data.head.idx.kspace_encode_step_1) == iY(idxY);
                        if(nnz(lIdx) > 1)
                            lIdx(find(lIdx,nnz(lIdx)-1,'last')) = false;
                        end
                        dataLine = data.data{lIdx};
                        
                        for idxCha = 1:measPara.dim(5)
                            kSpace{idxSli,idxCha,idxRep,idxAvg}(idxY,:,idxPha) = dataLine(2*iFreq*(idxCha-1)+1:2:2*iFreq*idxCha).'+ 1i * dataLine(2*iFreq*(idxCha-1)+2:2:2*iFreq*idxCha).';
                        end
                    end
                end
            end

        elseif(strcmp(measPara.dimension,'3D') || strcmp(measPara.dimension,'4D'))
            for idxZ = 1:measPara.dim(3)                               
                for idxPha = 1:measPara.dim(4)                        
                    for idxY = 1:measPara.dim(1)
                        lIdx = double(data.head.idx.repetition) == iRep(idxRep) & ...
                            double(data.head.idx.average) == iAvg(idxAvg) & ...
                            double(data.head.idx.kspace_encode_step_2) == iZ(idxZ) & ...
                            double(data.head.idx.phase) == iPha(idxPha) & ...
                            double(data.head.idx.kspace_encode_step_1) == iY(idxY);
                        if(nnz(lIdx) > 1)
                            lIdx(find(lIdx,nnz(lIdx)-1,'last')) = false;
                        end
                        dataLine = data.data{lIdx};
                        
                        for idxCha = 1:measPara.dim(5)
                            kSpace{1,idxCha,idxRep,idxAvg}(idxY,:,idxZ,idxPha) = dataLine(2*iFreq*(idxCha-1)+1:2:2*iFreq*idxCha).' + 1i * dataLine(2*iFreq*(idxCha-1)+2:2:2*iFreq*idxCha).';
                        end
                    end
                end
            end
        end
    end
end
end


