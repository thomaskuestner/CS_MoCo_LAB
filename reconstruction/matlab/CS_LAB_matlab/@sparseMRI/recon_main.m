function imageCha = recon_main(obj, iRep, iAvg)
% main reconstruction function
% prepare constraints and select algorithm depending on image dimensionality
%
% input:
% obj       CS reconstruction object (holding all parameters)
% iRep      current repetition
% iAvg      current average
%
% output:
% imageCha  individual reconstructed channel images
%
% (c) Thomas Kuestner 
% -------------------------------------------------------------------------

nCha = obj.measPara.dim(5);
nSlices = size(obj.kSpace,1);
nTime = obj.measPara.dim(4);
imageCha = cell(nSlices,nCha);

if(strcmp(obj.measPara.dimension,'2D'))
% lMask = cellfun(@(x) abs(x) > 0, all_KSpaces, 'UniformOutput', false);

    dispProgress('Slices', 0, nSlices);
    for iSli = 1:nSlices
        % prepare k-space 
        % => k_y - k_x - t - cha
        kSpaceL = cell2mat(shiftdim(obj.kSpace(iSli,:,iRep,iAvg),-2));
        if(strcmp(obj.measPara.precision,'double'))
            kSpaceL = double(kSpaceL);
        end
        % same undersampling for all coils and fully sampled in k_x direction
        % => reduce 4D problem (k_y - k_x - t - cha) to 2D problem (k_y - t)
        obj.fullMask = abs(kSpaceL(:,:,:,1)) > 0; % k_y - k_x - t
        if(obj.measPara.oversampling{2,1})
            kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - t - cha
            kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
            obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :);
            kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - t -cha
            kSpaceL = kSpaceL .* repmat(obj.fullMask,[1 1 1 size(kSpaceL,4)]);
            obj.measPara.dim(2) = size(kSpaceL,2);
        end
        [nPha, nFreq, nTime, nCha] = size(kSpaceL); 

        % prepare k-space => k_y - k_x - cha - t
        kSpaceL = permute(kSpaceL, [1 2 4 3]);
    
        % ensure that empty entries are zero
        mask = false(nPha,nFreq,nCha,nTime);
        for i=1:nTime
            mask(:,:,:,i) = repmat(obj.fullMask(:,:,i),[1,1,nCha]); 
        end
        kSpaceL(~mask) = 0;
        
        % set reconstruction parameter
        SMparam.TVWeight = obj.lambdaTV;
        SMparam.xfmWeight = obj.lambda;
        SMparam.Itnlim = obj.iNINNER;
        SMparam.l1Smooth = obj.l1Smooth;
        SMparam.pNorm = obj.p;
        SMparam.lineSearchItnlim = obj.lineSearchItnlim;
        SMparam.lineSearchAlpha = obj.lineSearchAlpha;
        SMparam.lineSearchBeta = obj.lineSearchBeta;
        SMparam.lineSearchT0 = obj.lineSearchT0; 
        
        if(strcmp(obj.trafo.trafoType,'fft'))
            SMParam.XFM = 1; 
        else % wavelet_lab
            % find the closest diadic size for the images
            ssx = 2^ceil(log2(size(kSpaceL,1))); 
            ssy = 2^ceil(log2(size(kSpaceL,2)));
            ss = max(ssx, ssy);
            SMParam.XFM = Wavelet(obj.trafo.waveletFilter, obj.trafo.waveletFilterSize, log2(ss/2^(obj.trafo.waveletStages)));
        end

        image = zeros(nPha,nFreq,nCha,nTime);
        dispProgress('Time', 0, nTime);
        fprintf('sparseMRI reconstruction\n');
        for i=1:nTime
            dispProgrss('Channel',0,nCha);
            for iCha=1:nCha
                if(obj.window.windowingOn)
                    window = windowND(obj.window.type,[size(kSpaceL,1), size(kSpaceL,2)],obj.window.windowOpt{1},obj.window.windowOpt{2});
                    phaseMap = exp(i*angle((ifft2c(kSpaceL(:,:,iCha,i).*window))));
                    clear 'window';
                else
                    phaseMap = 1;
                end
                param.FT = p2DFT(mask(:,:,iCha,i), [size(kSpaceL,1), size(kSpaceL,2)], phaseMap, 2); 
                
                res = zeros(size(kSpaceL,1),size(kSpaceL,2));
                for iO=1:iNOUTER
                    res = fnlCg(res, SMParam);
                end
                image(:,:,iCha,i) = SMParam.XFM'*res;
                dispProgress('Channel',iCha/nCha);
            end
            dispProgress('Channel','Close');
            dispProgress('Time',i/nTime);
        end
        dispProgress('Time','Close');
        
        imageCha(iSli,:) = shiftdim(mat2cell(permute(image,[1 2 4 3]),nPha,nFreq,nTime,ones(1,nCha)),2);

        dispProgress('Slices', iSli/nSlices);
    end
    dispProgress('Slices', 'Close');

elseif(strcmp(obj.measPara.dimension,'3D'))
        
    % prepare k-space 
    % => k_y - k_z - k_x - cha
    kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-2));
    
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    % same undersampling for all coils and fully sampled in k_x direction
    obj.fullMask = abs(kSpaceL) > 0; % k_y - k_x - k_z - cha
    if(obj.measPara.oversampling{2,1})
        kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - k_z - cha
        kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1},:, :); % anti-aliasing
        obj.fullMask = obj.fullMask(:,obj.measPara.oversampling{1,1},:, :);
        kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_z - k_x - cha
        kSpaceL = kSpaceL .* obj.fullMask;
        obj.measPara.dim(2) = size(kSpaceL,2);
    end
    
    kSpaceL = permute(kSpaceL,[1 3 2 4]);
    obj.fullMask = permute(obj.fullMask, [1 3 2 4]);
    [nPha, nZ, nFreq, nCha] = size(kSpaceL);   

    fprintf('Performing calibration\n');
    % perform calibration
    % get calibration data from k-space
    if(isempty(obj.calibSize))
        % automatically determine calibration region
        obj.calibSize = obj.calibrationSize(obj.fullMask(:,:,:,1)); % k_y - k_z - k_x
%         obj.kCalib = crop(kSpaceL, [obj.calibSize, nCha]); % k_y - k_z - k_x - cha
    else
        % fixed calibration region
        plonk = [nPha, nZ, nFreq];
        idx = cell(1,length(plonk));
        for iInner=1:length(plonk)
            if(mod(obj.calibSize(iInner),2) == 0)
                idx{iInner} = floor(plonk(iInner)/2)+1+ceil(-obj.calibSize(iInner)/2):floor(plonk(iInner)/2)+ceil(obj.calibSize(iInner)/2);
            else
                idx{iInner} = floor(plonk(iInner)/2)+ceil(-obj.calibSize(iInner)/2):floor(plonk(iInner)/2)+ceil(obj.calibSize(iInner)/2)-1;
            end

            tmp = [idx{iInner}(1) <= 0, idx{iInner}(end) > plonk(iInner)];
            if(any(tmp))
                if(all(tmp)), error('crop(): Both index out of bounds'); end;
                hShift = [abs(idx{iInner}(1)) + 1, idx{iInner}(end) - plonk(iInner)];
                op = {idx{iInner}(1), idx{iInner}(end); '+', '-'; '<', '>'; plonk(iInner) + 1, 0};
                eval(sprintf('if(op{1,~tmp} %s hShift(tmp) %s %d), idx{iInner} = idx{iInner} %s hShift(tmp); else idx{iInner} = idx{iInner}(idx{iInner} %s %d);end;',op{2,helper}, op{3,helper}, op{4,helper}, op{2,helper}, op{3,~helper}, op{4,~helper}));
            end
        end
%         obj.kCalib = kSpaceL(idx{1},idx{2},idx{3},:); % k_y - k_z - k_x - cha
    end
           
	imageCha = obj.CSNLCG(kSpaceL);
        
            
elseif(strcmp(obj.measPara.dimension,'4D'))
    % TODO

end
end