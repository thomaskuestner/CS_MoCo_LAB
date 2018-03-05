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

        
        %     obj.fullMask = lMask{1}(:, obj.obj.measPara.oversampling{1}, :); % k_y - k_x - t
        %     mask = abs(sum(permute(fullMask,[3 1 2]),3)) > 0;

        if(iSli == 1)
            if(isempty(obj.calibSize))
                obj.calibSize = cell(1,nTime);
                reempty = true;
            else
                % adjust calibSize due to anti-aliasing
                obj.calibSize(2) = nFreq;
                reempty = false;
            end
        else
            if(reempty)
                obj.calibSize = cell(1,nTime);
            end
        end
        obj.kCalib = cell(1,nTime);
        obj.kernel = obj.kCalib;
        obj.kernelImg = obj.kCalib;
        allGOPs = obj.kCalib;
        obj.A = obj.kCalib;
    
        % convolution with GRAPPA kernel in k-space -> multiplication in image
        % domain
        obj.method = 'fft';

        % prepare k-space => k_y - k_x - cha - t
        kSpaceL = permute(kSpaceL, [1 2 4 3]);
    
        fprintf('Performing calibration\n');
        for i=1:nTime

            % perform calibration
            % get calibration data from k-space
            if(iscell(obj.calibSize))
                % automatically determine calibration region
                obj.calibSize{i} = obj.calibrationSize(obj.fullMask(:,:,i)); % k_y - k_x (cell: t)
                obj.kCalib{i} = crop(kSpaceL(:,:,:,i), [obj.calibSize{i}, nCha]); % k_y - k_x - cha (cell: t)
            else
                % fixed calibration region
                plonk = [nPha, nFreq];
                idx = cell(1,length(plonk));
                for iInner=1:length(plonk)
                    idx{iInner} = floor(plonk(iInner)/2)+1+ceil(-obj.calibSize(iInner)/2):floor(plonk(iInner)/2)+ceil(obj.calibSize(iInner)/2);
                end
                obj.kCalib{i} = kSpaceL(idx{1},idx{2},:,i);
            end

%             % deprecated
%             obj.kernel{i} = zeros([obj.kernelSize, nCha, nCha]); % k_y - k_x - cha - cha (cell: t)
% 
%             [AtA,obj.A{i}] = obj.corrMatrix(i);
%             for n=1:nCha
%                 obj.kernel{i}(:,:,:,n) = obj.calibrate(AtA,nCha,n);
%             end
            
            obj.kernel = calibSPIRiT(obj.kCalib{i}, obj.kernelSize, nCha, obj.calibTyk);
            allGOPs{i} = SPIRiT(obj.kernel, 'fft',[nPha,nFreq]);

        end

        % ensure that empty entries are zero
        mask = false(nPha,nFreq,nCha,nTime);
        for i=1:nTime
            mask(:,:,:,i) = repmat(obj.fullMask(:,:,i),[1,1,nCha]); 
        end
        kSpaceL(~mask) = 0;

        image = zeros(nPha,nFreq,nCha,nTime);
        dispProgress('Time', 0, nTime);
        fprintf('SPIRiT reconstruction\n');
        for i=1:nTime
%             image(:,:,:,i) = sum(abs(ifft2c(obj.pocsSPIRiT(kSpaceL(:,:,:,i),allGOPs{i},kSpaceL(:,:,:,i),0))).^2,3); % -> y - x - t
            if(strcmp(obj.solver,'POCS'))
                image(:,:,:,i) = ifftnshift(obj.pocsSPIRiT(kSpaceL(:,:,:,i),allGOPs{i},kSpaceL(:,:,:,i),0),1:2); % => y-x-cha-t
            else
                image(:,:,:,i) = ifftnshift(obj.cgSPIRiT(kSpaceL(:,:,:,i),allGOPs{i},obj.iNINNER,obj.reconTyk,kSpaceL(:,:,:,i)),1:2);
            end
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
    kSpaceL = permute(kSpaceL,[1 3 2 4]);
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    % same undersampling for all coils and fully sampled in k_x direction
    obj.fullMask = abs(kSpaceL(:,:,:,1)) > 0; % k_y - k_z - k_x
    if(obj.measPara.oversampling{2,1})
        kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 3),[],3),3); % -> k_y - k_z - x - cha
        kSpaceL = kSpaceL(:, :, obj.measPara.oversampling{1,1}, :); % anti-aliasing
        obj.fullMask = obj.fullMask(:,:, obj.measPara.oversampling{1,1}, :);
        kSpaceL = fftshift(fft(ifftshift(kSpaceL,3),[],3),3); % -> k_y - k_z - k_x - cha
        kSpaceL = kSpaceL .* repmat(obj.fullMask, [1 1 1 size(kSpaceL,4)]);
        obj.measPara.dim(2) = size(kSpaceL,3);
    end
    
    [nPha, nZ, nFreq, nCha] = size(kSpaceL); 
    
    if(~isempty(obj.calibSize))
        % adjust calibSize due to anti-aliasing
        obj.calibSize(3) = nFreq;
    end
    
    % create a 3D kernel fomr the calibration region
    fprintf('Performing calibration\n');
    if(isempty(obj.calibSize))
        % automatically determine calibration region
        obj.calibSize = obj.calibrationSize(obj.fullMask(:,:,:,1)); % k_y - k_z - k_x
        if(length(obj.calibSize) == 2 && ndims(kSpaceL) > 2), obj.calibSize = [obj.calibSize(1), obj.calibSize(2), 1]; end;
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
                eval(sprintf('if(op{1,~helper} %s hShift(helper) %s %d), idx{iInner} = idx{iInner} %s hShift(helper); else idx{iInner} = idx{iInner}(idx{iInner} %s %d);end;',op{2,helper}, op{3,helper}, op{4,helper}, op{2,helper}, op{3,~helper}, op{4,~helper}));
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
    end
    
%     kernel = zeros([obj.kernelSize, nCha, nCha]); % k_y - k_z - k_x - cha - cha
%     kernelImg = zeros([nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha, nCha]);
        
    kernel3D = zeros([obj.kernelSize, nCha, nCha]); % k_y - k_z - k_x - cha - cha
    hybridKernel = zeros([obj.kernelSize(1), obj.kernelSize(2), nFreq+obj.kernelSize(3)-1, nCha, nCha]);
%     kernelImg3D = zeros([nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha, nCha]);
        
    if(nFreq == 1)
        [AtA, obj.A] = obj.corrMatrix();
    else
        [AtA,obj.A] = obj.corrMatrix3D();
    end
    for n=1:nCha
        if(nFreq == 1)
            kernel3D(:,:,:,:,n) = obj.calibrate(AtA,nCha,n);
        else
            kernel3D(:,:,:,:,n) = obj.calibrate3D(AtA,nCha,n);
        end

        hybridKernel(:,:,:,:,n) = ifftnshift(zpad(kernel3D(:,:,end:-1:1,:,n), obj.kernelSize(1), obj.kernelSize(2), nFreq+obj.kernelSize(3)-1, nCha),3);

        % convert correlation kernel in k-space to convolution kernel in
        % image space
        % zpad for circular conv = linear conv
%         kernelImg3D(:,:,:,:,n) = ifftnshift(zpad(kernel3D(end:-1:1,end:-1:1,end:-1:1,:,n), nPha+obj.kernelSize(1)-1, nZ+obj.kernelSize(2)-1, nFreq+obj.kernelSize(3)-1, nCha)); % y - z - x - cha - cha
    end
        
    % compute inverse fourier transform along readout direction
    kSpaceL = ifftnshift(kSpaceL,3); % k_y - k_z - x - cha
    
    % convolution with GRAPPA kernel in k-space -> multiplication in image
    % domain
    obj.method = 'fft';

    image = zeros(nPha,nFreq,nZ,nCha);
%     AtA_sub = zeros(prod(obj.kernelSize)*nCha); % subtract this matrix from AtA_mask, due to accumulation of previous samples (A'*A multiplication)
%     A_sub = zeros((size(obj.kCalib,1)-obj.kernelSize(1)+1) * (size(obj.kCalib,2)-obj.kernelSize(2)+1) * (size(obj.kCalib,3)-obj.kernelSize(3)+1),prod(obj.kernelSize)*nCha);

    fprintf('SPIRiT reconstruction\n');
    dispProgress('Frequency',0,nFreq);
    for iFreq=1:nFreq
        
        % extract 2D kernel from complete 3D kernel
        obj.kernel = squeeze(hybridKernel(:,:,iFreq,:,:));
  %         maskHelper = zeros(size(obj.kCalib));
%         maskHelper(:,:,iFreq,:) = 1;
%         obj.kCalib = maskHelper;
%         [~, A_mask] = obj.corrMatrix3D();
%         AtA_mask = AtA_mask - AtA_sub;
%         AtA_sub = AtA_sub + AtA_mask;
%         AtA_mask = logical(AtA_mask);
%         AtA_2D = reshape(AtA(AtA_mask),obj.kernelSize(1)*obj.kernelSize(2)*nCha,obj.kernelSize(1)*obj.kernelSize(2)*nCha);
%         A_mask = A_mask - A_sub;
%         A_sub = A_sub + A_mask;
%         A_mask = logical(A_mask);
%         A_2D = reshape(obj.A(A_mask),(size(obj.kCalib,1)-obj.kernelSize(1)+1) * (size(obj.kCalib,2)-obj.kernelSize(2)+1),obj.kernelSize(1)*obj.kernelSize(2)*nCha);
%         AtA_2D = A_2D' * A_2D;
        
%         obj.kernelImg = squeeze(kernelImg3D(:,:,iFreq,:));
%         obj.kernel = zeros([obj.kernelSize(1), obj.kernelSize(2), nCha, nCha]); % k_y - k_z - cha - cha
%         for n=1:nCha
% %             obj.kernel(:,:,:,n) = obj.calibrate(AtA_2D,nCha,n);
%         end
        
%         % perform calibration
%         % get calibration data from k-space
%         if(isempty(obj.calibSize))
%             % automatically determine calibration region
%             obj.calibSize = obj.calibrationSize(obj.fullMask(:,:,iFreq)); % k_y - k_z
%             obj.kCalib = crop(squeeze(kSpaceL(:,:,iFreq,:)), [obj.calibSize, nCha]); % k_y - k_z - cha
%         else
%             % fixed calibration region
%             plonk = [nPha, nZ];
%             idx = cell(1,length(plonk));
%             for iInner=1:length(plonk)
%                 if(mod(obj.calibSize(iInner),2) == 0)
%                     idx{iInner} = floor(plonk(iInner)/2)+1+ceil(-obj.calibSize(iInner)/2):floor(plonk(iInner)/2)+ceil(obj.calibSize(iInner)/2);
%                 else
%                     idx{iInner} = floor(plonk(iInner)/2)+ceil(-obj.calibSize(iInner)/2):floor(plonk(iInner)/2)+ceil(obj.calibSize(iInner)/2)-1;
%                 end
% 
%                 tmp = [idx{iInner}(1) <= 0, idx{iInner}(end) > plonk(iInner)];
%                 if(any(tmp))
%                     if(all(tmp)), error('crop(): Both index out of bounds'); end;
%                     hShift = [abs(idx{iInner}(1)) + 1, idx{iInner}(end) - plonk(iInner)];
%                     op = {idx{iInner}(1), idx{iInner}(end); '+', '-'; '<', '>'; plonk(iInner) + 1, 0};
%                     eval(sprintf('if(op{1,~helper} %s hShift(helper) %s %d), idx{iInner} = idx{iInner} %s hShift(helper); else idx{iInner} = idx{iInner}(idx{iInner} %s %d);end;',op{2,helper}, op{3,helper}, op{4,helper}, op{2,helper}, op{3,~helper}, op{4,~helper}));
%                 end
%             end
%             obj.kCalib = squeeze(kSpaceL(idx{1},idx{2},iFreq,:));
%         end
% 
%         obj.kernel = zeros([obj.kernelSize(1), obj.kernelSize(2), nCha, nCha]); % k_y - k_z - cha - cha
% 
%         [AtA,obj.A] = obj.corrMatrix();
%         for n=1:nCha
%             obj.kernel(:,:,:,n) = obj.calibrate(AtA,nCha,n);
%         end

        allGOPs = SPIRiT(obj.kernel, 'fft',[nPha,nZ]);
        
        if(strcmp(obj.solver,'POCS'))
            image(:,iFreq,:,:) = permute(obj.pocsSPIRiT(squeeze(kSpaceL(:,:,iFreq,:)),allGOPs,squeeze(kSpaceL(:,:,iFreq,:)),0),[1 4 2 3]);
        else
            image(:,iFreq,:,:) = permute(obj.cgSPIRiT(squeeze(kSpaceL(:,:,iFreq,:)),allGOPs,obj.iNINNER,obj.reconTyk,squeeze(kSpaceL(:,:,iFreq,:))),[1 4 2 3]);
        end
        
        dispProgress('Frequency', iFreq/nFreq);
    end
    dispProgress('Frequency','Close');
    image = ifftnshift(image,[1,3]);
    imageCha = shiftdim(mat2cell(image,nPha,nFreq,nZ,ones(1,nCha)),2);
            
elseif(strcmp(obj.measPara.dimension,'4D'))
    % TODO

end
end