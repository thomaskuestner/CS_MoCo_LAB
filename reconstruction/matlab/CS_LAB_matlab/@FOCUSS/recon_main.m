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
% nTime = obj.measPara.dim(4);
imageCha = cell(nSlices,nCha);

if(strcmp(obj.measPara.dimension,'2D'))
    % 2D case: x-y-t or x-y (sparse dimensions: y-t)
    
    dispProgress('Slices', 0, nSlices);
    for iSli = 1:nSlices      
        % prepare k-space 
        % => k_y - k_x - t - cha
        kSpaceL = cell2mat(shiftdim(obj.kSpace(iSli,:,iRep,iAvg),-2));
        if((obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1 && nSlices == 1) || ...
                (iRep == size(obj.kSpace,3) && iAvg == size(obj.kSpace,4) && iSli == nSlices) || ...
                (obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1 && iSli == nSlices))
            obj.kSpace = [];
        end
        if(strcmp(obj.measPara.precision,'double'))
            kSpaceL = double(kSpaceL);
        end
        % extraction of sampling mask
        % !!! Must be the first operation, because fft and ifft operations,
        % destroy this information !!!
        % same undersampling for all coils and fully sampled in k_x direction
        obj.fullMask = abs(kSpaceL) > 0; % k_y - k_x - t - cha
        if(obj.measPara.oversampling{2,1})
            obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
            kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - t - cha
            kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
            kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - t -cha
            kSpaceL = kSpaceL .* obj.fullMask; % just take acquired points
            obj.measPara.dim(2) = size(kSpaceL,2);
        end
        
        [nPha, nFreq, nTime, nCha] = size(kSpaceL);

        if(iSli == 1)
            if(isempty(obj.calibSize))
                obj.calibSize = cell(1,nTime);
                reempty = true;
            else
                % adjust calibSize due to anti-aliasing
                if(~iscell(obj.calibSize))
                    obj.calibSize(2) = nFreq;
                    reempty = false;
                else % due to several averages
                    obj.calibSize = cell(1,nTime); 
                    reempty = true;
                end
            end
        else
            if(reempty)
                obj.calibSize = cell(1,nTime);
            end
        end
        obj.kCalib = cell(1,nTime);
        kernel = obj.kCalib;
        obj.kernelImg = obj.kCalib;

        % prepare k-space and mask => k_y - k_x - cha - t
        kSpaceL = permute(kSpaceL, [1 2 4 3]);
        obj.fullMask = permute(obj.fullMask, [1 2 4 3]);

        fprintf('Performing calibration\n');
        for i=1:nTime
            % perform calibration
            % get calibration data from k-space
            if(iscell(obj.calibSize))
                % automatically determine calibration region
                obj.calibSize{i} = obj.calibrationSize(obj.fullMask(:,:,1,i)); % k_y - k_x (cell: t)
                obj.kCalib{i} = crop(kSpaceL(:,:,:,i), [obj.calibSize{i}, nCha]); % k_y - k_x - cha (cell: t)
            else
                % fixed calibration region
                plonk = [nPha, nFreq];
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
                obj.kCalib{i} = kSpaceL(idx{1},idx{2},:,i);
            end

            % windowing for cut-out region
            if(obj.window.windowingOn)
                window = windowND(obj.window.type,[size(obj.kCalib{i},1), size(obj.kCalib{i},2)],obj.window.windowOpt{1},obj.window.windowOpt{2});
                obj.kCalib{i} = obj.kCalib{i} .* repmat(window,[1 1 nCha]);
                clear 'window';
            end

            if(obj.lambdaCalib > 0)
                kernel{i} = zeros([obj.kernelSize, nCha, nCha],obj.measPara.precision); % k_y - k_x - cha - cha (cell: t)
                if(obj.flagZeropadding)
                    obj.kernelImg{i} = zeros([nPha+obj.kernelSize(1)-1, nFreq+obj.kernelSize(2)-1, nCha, nCha],obj.measPara.precision); % k_y - k_x - cha - cha (cell: t)
                else
                    obj.kernelImg{i} = zeros([nPha, nFreq, nCha, nCha],obj.measPara.precision); % k_y - k_x - cha - cha (cell: t)
                end

                [AtA,~] = obj.corrMatrix2D(i);
                if(i == 1), dispProgress('GRAPPA calibration',0,nCha*nTime); end;
                for n=1:nCha
                    kernel{i}(:,:,:,n) = obj.calibrate2D(AtA,nCha,n);

                    % convert correlation kernel in k-space to convolution kernel in
                    % image space
                    % zpad for circular conv = linear conv
                    if(obj.flagZeropadding)
                        obj.kernelImg{i}(:,:,:,n) = (ifftnshift(zpad(kernel{i}(end:-1:1,end:-1:1,:,n), nPha+obj.kernelSize(1)-1, nFreq+obj.kernelSize(2)-1, nCha),1:2)); % y - x - cha - cha (cell: t)
                    else
                        obj.kernelImg{i}(:,:,:,n) = (ifftnshift(zpad(kernel{i}(end:-1:1,end:-1:1,:,n), nPha, nFreq, nCha),1:2)); % y - x - cha - cha (cell: t)
                    end
                    dispProgress('GRAPPA calibration',((i-1)* nCha + n)/(nCha*nTime));
                end
                if(i == nTime), dispProgress('GRAPPA calibration','Close'); end;
            end
        end
        
        % remove all unnecessary variables
        obj.kCalib = [];
        clear 'kernel' 'window' 'AtA';

        % k-t FOCUSS reconstruction
        fprintf('FOCUSS reconstruction\n');
        if(nTime == 1)
            imageCha(iSli,:) = obj.ktFOCUSS(kSpaceL);
        else
            imageCha(iSli,:) = obj.ktFOCUSS2Dt(kSpaceL);
        end

        dispProgress('Slices', iSli/nSlices);
    end
    dispProgress('Slices', 'Close');
    
elseif(strcmp(obj.measPara.dimension,'3D'))
    % 3D case: x-y-z (sparse dimensions: y-z)
           
    % prepare k-space 
    % => k_y - k_x - k_z - cha
    kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-2));
    if((obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1) || (iRep == size(obj.kSpace,3) && iAvg == size(obj.kSpace,4)))
        obj.kSpace = [];
    end
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    % extraction of sampling mask
    % !!! Must be the first operation, because fft and ifft operations,
    % destroy this information !!!
    % same subsampling for all channels
    obj.fullMask = abs(kSpaceL) > 0; % -> k_y-k_x-k_z-cha
    
    if(length(unique(sum(obj.fullMask,2))) > 2) % ensure same subsampling along readout direction
        val = max(squeeze(sum(sum(sum(obj.fullMask,1),2),3)))/size(obj.fullMask,2);
        i = 1;
        while(sum(sum(obj.fullMask(:,i,:,1))) ~= val && i <= size(obj.fullMask,2))
            i = i+1;
        end
        if(i > size(obj.fullMask,2))
            error('recon_main(): Incomplete undersampling mask!');
        end
        obj.fullMask = repmat(obj.fullMask(:,i,:,:), [1 size(obj.fullMask,2) 1 1]);
    end
        
    if(obj.measPara.oversampling{2,1})
        obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
        kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - k_z - cha
        kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
        kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - k_z - cha
        kSpaceL = kSpaceL .* obj.fullMask; % just take acquired points
        obj.measPara.dim(2) = size(kSpaceL,2);
    end
    
    [nPha, nFreq, nZ, nCha] = size(kSpaceL); 
    
    if(~isempty(obj.calibSize))
        % adjust calibSize due to anti-aliasing
        obj.calibSize(3) = nFreq;
    end
    
    % prepare k-space and mask => k_y - k_z - k_x - cha
    kSpaceL = permute(kSpaceL, [1 3 2 4]);
    obj.fullMask = permute(obj.fullMask, [1 3 2 4]);
    
    fprintf('Performing calibration\n');
    % perform calibration
    % get calibration data from k-space
    if(isempty(obj.calibSize))
        % automatically determine calibration region
        obj.calibSize = obj.calibrationSize(obj.fullMask(:,:,:,1)); % k_y - k_z - k_x
        if(length(obj.calibSize) == 2 && ndims(kSpaceL) > 2), obj.calibSize = [obj.calibSize, 1]; end;
        obj.kCalib = crop(kSpaceL, [obj.calibSize, nCha]); % k_y - k_z - k_x - cha
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
        obj.kCalib = kSpaceL(idx{1},idx{2},idx{3},:); % k_y - k_z - k_x - cha
    end

    % windowing for cut-out region
%     [~,w1] = window2D(prop.window.type,size(obj.kCalib,1),prop.window.windowOpt{1},prop.window.windowOpt{2});
%     [~,w2] = window2D(prop.window.type,size(obj.kCalib,2),prop.window.windowOpt{1},prop.window.windowOpt{2});
%     window =  w1(:) * w2(:).';
    if(obj.window.windowingOn)
        window = windowND(obj.window.type,[size(obj.kCalib,1), size(obj.kCalib,2)],obj.window.windowOpt{1},obj.window.windowOpt{2});
        obj.kCalib = obj.kCalib .* repmat(window,[1 1 size(obj.kCalib,3) nCha]);
        clear 'window';
    end

    if(obj.lambdaCalib > 0)
        kernel = zeros([obj.kernelSize, nCha, nCha],obj.measPara.precision); % k_y - k_z - k_x - cha - cha
        if(obj.flagZeropadding)
            obj.kernelImg = zeros([nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha, nCha],obj.measPara.precision);
        else
            obj.kernelImg = zeros([nPha, nZ, nFreq, nCha, nCha],obj.measPara.precision);
        end

        [AtA,~] = obj.corrMatrix3D();
        dispProgress('GRAPPA calibration',0,nCha);
        for n=1:nCha
            kernel(:,:,:,:,n) = obj.calibrate3D(AtA,nCha,n);

            % convert correlation kernel in k-space (conv) to convolution kernel in
            % image space (mult)
            % zpad for circular conv = linear conv
            if(obj.flagZeropadding)
                obj.kernelImg(:,:,:,:,n) = ifftnshift(zpad(kernel(end:-1:1,end:-1:1,end:-1:1,:,n), nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha),1:3); % y - z - x - cha - cha
            else
                obj.kernelImg(:,:,:,:,n) = ifftnshift(zpad(kernel(end:-1:1,end:-1:1,end:-1:1,:,n), nPha, nZ, nFreq, nCha),1:3); % y - z - x - cha - cha
            end
            dispProgress('GRAPPA calibration', n/nCha);
        end
        dispProgress('GRAPPA calibration', 'Close');
    end
    
    % remove all unnecessary variables
    obj.kCalib = [];
    clear 'kSpace' 'kernel' 'window' 'AtA';
    
    % k-t FOCUSS reconstruction
    fprintf('FOCUSS reconstruction\n');
    imageCha = obj.ktFOCUSS3D(kSpaceL);
    
    
elseif(strcmp(obj.measPara.dimension,'4D'))
    % 4D case: x-y-z-t (3D spatial and temporal, sparse dimensions: y-z-t)   
    
    % prepare k-space 
    % k-space: k_y - k_x - k_z - t - cha
    kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-3));
    if((obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1) || (iRep == size(obj.kSpace,3) && iAvg == size(obj.kSpace,4)))
        obj.kSpace = [];
    end
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    % same undersampling for all coils and fully sampled in k_x direction
    obj.fullMask = abs(kSpaceL) > 0; % k_y - k_z - k_x - t - cha
    if(obj.measPara.oversampling{2,1})
        obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :, :); % anti-aliasing
        kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - k_z - t - cha
        kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :, :); % anti-aliasing
        kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - k_z - t - cha
        kSpaceL = kSpaceL .* obj.fullMask; % just take acquired points
        obj.measPara.dim(2) = size(kSpaceL,2);
    end
    
    [nPha, nFreq, nZ, nTime, nCha] = size(kSpaceL);    
    
    if(isempty(obj.calibSize))
        lCalibauto = true;
        obj.calibSize = cell(1,nTime);
    else
        lCalibauto = false; % fixed calibration
    end
    obj.kCalib = cell(1,nTime);
    obj.kernelImg = obj.kCalib;

    
    % prepare k-space and mask => k_y - k_z - k_x - cha - t
    kSpaceL = permute(kSpaceL, [1 3 2 5 4]);
    obj.fullMask = permute(obj.fullMask, [1 3 2 5 4]);
    
    fprintf('Performing calibration');
    for i=1:nTime
        fprintf('.');
        
        % perform calibration
        % get calibration data from k-space
        if(lCalibauto)
            % automatically determine calibration region
            obj.calibSize{i} = obj.calibrationSize(obj.fullMask(:,:,:,1,i)); % k_y - k_z - k_x (cell: t)
            if(i >= 2) % normally first gate (end-exp.) is nearly full
                if(any(obj.calibSize{i} < 0.5 .* obj.calibSize{1}) || (obj.lambdaCalib > 0 && any(obj.calibSize{i} < obj. kernelSize)))
                    obj.calibSize{i} = obj.calibrationSize(obj.fullMask(:,:,:,1,i),[0.1 0.1 -1]); % allow some zero values inside the calibration region
                end
            end
            if((obj.lambdaCalib > 0 && any(obj.calibSize{i} < obj. kernelSize))) % problem still persists -> increase calibSize to minimal kernelSize
                obj.calibSize{i}(obj.calibSize{i} < obj.kernelSize) = obj.kernelSize(obj.calibSize{i} < obj.kernelSize);
            end
            if(obj.lambdaCalib > 0), obj.kCalib{i} = crop(kSpaceL(:,:,:,:,i), [obj.calibSize{i}, nCha]); end; % k_y - k_z - k_x - cha (cell: t)
        else
            % fixed calibration region
            % update calibration size (starting from fixed region)
            if(i==1)
                fixedCalibSize = obj.calibSize;
                obj.calibSize = cell(1,nTime);
            end            
            obj.calibSize{i} = obj.calibrationSize(obj.fullMask(:,:,:,1,i),[0.075 0.1 -1], fixedCalibSize);
            if(obj.lambdaCalib > 0), obj.kCalib{i} = crop(kSpaceL(:,:,:,:,i), [obj.calibSize, nCha]); end;
        end
        
        % windowing for cut-out region
        if(obj.window.windowingOn)
            window = windowND(obj.window.type,[size(obj.kCalib{i},1), size(obj.kCalib{i},2)],obj.window.windowOpt{1},obj.window.windowOpt{2});
            obj.kCalib{i} = obj.kCalib{i} .* repmat(window,[1 1 size(obj.kCalib{i},3) nCha]);
            clear 'window';
        end
        
        if(obj.lambdaCalib > 0)
            kernel = zeros([obj.kernelSize, nCha, nCha],obj.measPara.precision); % k_y - k_z - k_x - cha - cha (cell: t)
            if(obj.flagZeropadding)
                obj.kernelImg{i} = zeros([1,nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha, nCha],obj.measPara.precision); % k_y - k_z - k_x - cha - cha (cell: t)
            else
                obj.kernelImg{i} = zeros([1,nPha, nZ, nFreq, nCha, nCha],obj.measPara.precision); % k_y - k_z - k_x - cha - cha (cell: t)
            end

            [AtA,~] = obj.corrMatrix4D(i);
            dispProgress('GRAPPA calibration',0,nCha);
            for n=1:nCha
                kernel(:,:,:,:,n) = obj.calibrate4D(AtA,nCha,n);

                % convert correlation kernel in k-space to convolution kernel in
                % image space
                % zpad for circular conv = linear conv
                if(obj.flagZeropadding)
                    obj.kernelImg{i}(1,:,:,:,:,n) = ifftnshift(zpad(kernel(end:-1:1,end:-1:1,end:-1:1,:,n), nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha),1:3); % y - z - x - cha - cha (cell: t)
                else
                    obj.kernelImg{i}(1,:,:,:,:,n) = ifftnshift(zpad(kernel(end:-1:1,end:-1:1,end:-1:1,:,n), nPha, nZ, nFreq, nCha),1:3); % y - z - x - cha - cha (cell: t)
                end
                dispProgress('GRAPPA calibration',n/nCha);
            end
            dispProgress('GRAPPA calibration','Close');
        end
    end
    fprintf('\n');
    obj.kernelImg = cell2mat(shiftdim(obj.kernelImg,1));
    
    % remove all unnecessary variables
    obj.kCalib = [];
    clear 'kSpace' 'kernel' 'window' 'AtA';
    
    kSpaceL = permute(kSpaceL, [5 1 2 3 4]); % => t - k_y - k_z - k_x - cha
    if(any(obj.trafo.fftBA(1,:))) % do transformations here to avoid lazy-copy in matlab when handing over kSpace variable
        for iBA=find(obj.trafo.fftBA(1,:))
            kSpaceL = ifftnshift(kSpaceL,iBA);
        end
    end
    
    % k-t FOCUSS reconstruction
    fprintf('FOCUSS reconstruction\n');
    if(obj.lambdaMC == 0)
        imageCha = obj.ktFOCUSS4D(kSpaceL);
    else
        imageCha = obj.ktFOCUSS4DMotion(kSpaceL);
    end

elseif(strcmp(obj.measPara.dimension,'5D'))
    % 5D case: x-y-z-t (3D spatial, temporal and gating domain, sparse dimensions: y-z-t-g)   
    
    % prepare k-space 
    % k-space: k_y - k_x - k_z - t - cha
    if(iscell(obj.kSpace))
        kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-4));
    else
        kSpaceL = obj.kSpace;
        imageCha = cell(1,size(kSpaceL,6));
    end
    if(obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1)
        obj.kSpace = [];
    end
    % same undersampling for all coils and fully sampled in k_x direction
    obj.fullMask = abs(kSpaceL) > 0; % k_y - k_z - k_x - t - g - cha
    if(obj.measPara.oversampling{2,1})
        obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :, :, :); % anti-aliasing

        kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - k_z - t - cha
        kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :, :, :); % anti-aliasing
        kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - k_z - t - cha
        kSpaceL = kSpaceL .* obj.fullMask; % just take acquired points
        obj.measPara.dim(2) = size(kSpaceL,2);
    end
    
    [nPha, nFreq, nZ, nTime, nSmp, nCha] = size(kSpaceL);    
    
    if(isempty(obj.calibSize))
        obj.calibSize = cell(1,nTime);
    end
    obj.kCalib = cell(1,nTime);
    obj.kernelImg = obj.kCalib;
    
    % prepare k-space and mask => k_y - k_z - k_x - cha - t - g
    kSpaceL = permute(kSpaceL, [1 3 2 6 4 5]);
    obj.fullMask = permute(obj.fullMask, [1 3 2 6 4 5]);
    
    fprintf('\nsize(img): ');
    fprintf('%d ', size(kSpaceL));
    fprintf('\n');
    
    fprintf('Performing calibration');
    for i=1:nTime
        fprintf('.');
        
        % perform calibration
        % get calibration data from k-space
        if(iscell(obj.calibSize))
            % automatically determine calibration region
            obj.calibSize{i} = obj.calibrationSize(obj.fullMask(:,:,:,1,i,1)); % k_y - k_z - k_x (cell: t)
            tmpSize = round(0.15 .* [nPha, nZ, nFreq]);
            if(any(obj.calibSize{i} < tmpSize))
                lIdx = obj.calibSize{i} < tmpSize;
                obj.calibSize{i}(lIdx) = subsref(tmpSize,struct('type','()','subs',{{lIdx}}));
            end
            obj.kCalib{i} = crop(kSpaceL(:,:,:,:,i,1), [obj.calibSize{i}, nCha]); % k_y - k_z - k_x - cha (cell: t)
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
            obj.kCalib{i} = kSpaceL(idx{1},idx{2},idx{3},:,i,:);
        end
        
        % windowing for cut-out region
        if(obj.window.windowingOn)
            window = windowND(obj.window.type,[size(obj.kCalib{i},1), size(obj.kCalib{i},2)],obj.window.windowOpt{1},obj.window.windowOpt{2});
            obj.kCalib{i} = obj.kCalib{i} .* repmat(window,[1 1 size(obj.kCalib{i},3) nCha]);
        end
        
        if(obj.lambdaCalib > 0)
            kernel = zeros([obj.kernelSize, nCha, nCha]); % k_y - k_z - k_x - cha - cha (cell: t)
            if(obj.flagZeropadding)
                obj.kernelImg{i} = zeros([nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha, nCha]); % k_y - k_z - k_x - cha - cha (cell: t)
            else
                obj.kernelImg{i} = zeros([nPha, nZ, nFreq, nCha, nCha]); % k_y - k_z - k_x - cha - cha (cell: t)
            end

            [AtA,~] = obj.corrMatrix4D(i);
            dispProgress('GRAPPA calibration',0,nCha);
            for n=1:nCha
                kernel(:,:,:,:,n) = obj.calibrate4D(AtA,nCha,n);

                % convert correlation kernel in k-space to convolution kernel in
                % image space
                % zpad for circular conv = linear conv
                if(obj.flagZeropadding)
                    obj.kernelImg{i}(:,:,:,:,n) = ifftnshift(zpad(kernel(end:-1:1,end:-1:1,end:-1:1,:,n), nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha),1:3); % y - z - x - cha - cha (cell: t)
                else
                    obj.kernelImg{i}(:,:,:,:,n) = ifftnshift(zpad(kernel(end:-1:1,end:-1:1,end:-1:1,:,n), nPha, nZ, nFreq, nCha),1:3); % y - z - x - cha - cha (cell: t)
                end
                dispProgress('GRAPPA calibration',n/nCha);
            end
            dispProgress('GRAPPA calibration','Close');
        end
    end
    fprintf('\n');
    
    % remove all unnecessary variables
    obj.kCalib = [];
    clear 'kSpace' 'kernel' 'window' 'AtA';
    
    kSpaceL = permute(kSpaceL, [5 1 2 3 6 4]); % => t - k_y - k_z - k_x - g - cha
    if(any(obj.trafo.fftBA(1,:))) % do transformations here to avoid lazy-copy in matlab when handing over kSpace variable
        for iBA=find(obj.trafo.fftBA(1,:))
            kSpaceL = ifftnshift(kSpaceL,iBA);
        end
    end
    
    % k-t FOCUSS reconstruction
    fprintf('FOCUSS reconstruction\n');
    imageCha = obj.ktFOCUSS5D(kSpaceL);
    
    
end

end