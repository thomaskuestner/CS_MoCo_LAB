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


elseif(strcmp(obj.measPara.dimension,'3D'))
    
    % prepare k-space 
    % => k_y - k_z - k_x - cha
    kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-2));
    kSpaceL = permute(kSpaceL,[1 3 2 4]);
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    % same undersampling for all coils and fully sampled in k_x direction
    obj.fullMask = abs(kSpaceL) > 0; % k_y - k_z - k_x
    if(obj.measPara.oversampling{2,1})
        kSpaceL = ifftnshift(kSpaceL, 3); % -> k_y - k_z - x - cha
        kSpaceL = kSpaceL(:, :, obj.measPara.oversampling{1,1}, :); % anti-aliasing
        obj.fullMask = obj.fullMask(:,:, obj.measPara.oversampling{1,1}, :);
        kSpaceL = fftnshift(kSpaceL,3); % -> k_y - k_z - k_x - cha
        kSpaceL = kSpaceL .* obj.fullMask;
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
        if(size(kSpaceL,3) == 1), obj.calibSize = [obj.calibSize, 1]; obj.kernelSize = obj.kernelSize([1 3 2]); end; % k_y - k_x sparse!
        if(~obj.flagToolbox)
            obj.kCalib = crop(kSpaceL, [obj.calibSize, nCha]); % k_y - k_z - k_x - cha
        end
    else
        if(~obj.flagToolbox)
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
    end
    
    if(~obj.flagToolbox)
        % windowing for cut-out region
    %     [~,w1] = window2D(prop.window.type,size(obj.kCalib,1),prop.window.windowOpt{1},prop.window.windowOpt{2});
    %     [~,w2] = window2D(prop.window.type,size(obj.kCalib,2),prop.window.windowOpt{1},prop.window.windowOpt{2});
    %     window =  w1(:) * w2(:).';
        if(obj.window.windowingOn)
            window = windowND(obj.window.type,[size(obj.kCalib,1), size(obj.kCalib,2)],obj.window.windowOpt{1},obj.window.windowOpt{2});
            obj.kCalib = obj.kCalib .* repmat(window,[1 1 size(obj.kCalib,3) nCha]);
        end

        % compute eigen-value maps
        [hybridKernel,S] = dat3Kernel(obj.kCalib,obj.kernelSize);
        hybridKernel = ifftnshift(zpad(hybridKernel(:,:,end:-1:1,:,:),[obj.kernelSize(1),obj.kernelSize(2),nFreq + obj.kernelSize(3)-1,nCha,size(hybridKernel,5)]),3);
        idx = max(find(S >= S(1)*obj.eigThresh_k));

        % compute inverse fourier transform along readout direction
        kSpaceL = ifftnshift(kSpaceL,3); % k_y - k_z - x - cha

        % convolution with GRAPPA kernel in k-space -> multiplication in image
        % domain
        obj.method = 'fft';

        % prepare wavelet transformation
        if(strcmp(obj.solver,'L1'))
            if(strcmp(obj.trafo.trafoType,'wavelet_lab'))
                ssx = 2^ceil(log2(nPha)); 
                ssy = 2^ceil(log2(nZ));
                ss = max(ssx, ssy);
                XOP = Wavelet(obj.trafo.waveletFilter,obj.trafo.waveletFilterSize,log2(ss/2^(obj.trafo.waveletStages)));
            end
        end

        image = zeros(nPha,nFreq,nZ,obj.n_maps);
    %     AtA_sub = zeros(prod(obj.kernelSize)*nCha); % subtract this matrix from AtA_mask, due to accumulation of previous samples (A'*A multiplication)
    %     A_sub = zeros((size(obj.kCalib,1)-obj.kernelSize(1)+1) * (size(obj.kCalib,2)-obj.kernelSize(2)+1) * (size(obj.kCalib,3)-obj.kernelSize(3)+1),prod(obj.kernelSize)*nCha);

        fprintf('ESPIRiT reconstruction\n');
        dispProgress('Frequency',0,nFreq);
        for iFreq=1:nFreq

            % extract 2D kernel from complete 3D kernel
            obj.kernel = squeeze(hybridKernel(:,:,iFreq,:,:));

            % crop kernels and compute eigen-value decomposition in image space to get maps
            [M,W] = kernelEig(obj.kernel(:,:,:,1:idx),[nPha,nZ]);

            % Compute Soft-SENSE ESPIRiT Maps 
            % crop sensitivity maps according to eigenvalues == 1
            % weigth with n maps of eigen-values

            maps = M(:,:,:,end-(obj.n_maps-1):end);

            % Weight the eigenvectors with soft-senses eigen-values
            weights = W(:,:,end-(obj.n_maps-1):end);
            weights = (weights - obj.eigThresh_im)./(1-obj.eigThresh_im).* (W(:,:,end-(obj.n_maps-1):end) > obj.eigThresh_im);
            weights = -cos(pi*weights)/2 + 1/2;

            % create an ESPIRiT operator
            ESP = ESPIRiT(maps,weights);

            if(strcmp(obj.solver,'L1'))
                FT = p2DFT(squeeze(obj.fullMask(:,:,iFreq,:)),[nPha,nZ,nCha]);
                if(strcmp(obj.trafo.trafoType,'fft'))
                    XOP = FT;
                end
                image(:,iFreq,:,:) = permute(obj.cgL1ESPIRiT(squeeze(kSpaceL(:,:,iFreq,:)), zeros([size(squeeze(kSpaceL(:,:,iFreq,1))),obj.n_maps]), FT, ESP, obj.iNINNER, XOP, obj.trafo.wavWeight, obj.splitWeight, obj.iNIterSplit),[1 4 2 3]); % => y-x-z-cha
            else
                [tmpK,tmp] = obj.cgESPIRiT(squeeze(kSpaceL(:,:,iFreq,:)), ESP, obj.iNINNER, obj.calibTyk, zeros(size(squeeze(kSpaceL(:,:,iFreq,:)))));
                image(:,iFreq,:,:) = permute(tmp,[1 4 2 3]);
            end
            dispProgress('Frequency', iFreq/nFreq);
        end
        dispProgress('Frequency','Close');
    %     image = ifftnshift(image,[1,3]);
        imageCha(:,1:obj.n_maps) = shiftdim(mat2cell(image,nPha,nFreq,nZ,ones(1,obj.n_maps)),2);
        [imageCha{:,obj.n_maps+1:end}] = deal(zeros(size(imageCha{1,1})));
    
    else
        fprintf('ESPIRiT reconstruction\n');
        if(ispc)
            eval(sprintf('ESP = bart(''ecalib -r %d:%d -k %d:%d -m %d'',kSpaceL);', obj.calibSize(1), obj.calibSize(2), obj.kernelSize(1), obj.kernelSize(2), obj.n_maps ));
            if(strcmp(obj.solver,'L1'))
                eval(sprintf('image = permute(squeeze(bart(''pics -l1 -r %f'', kSpaceL, ESP)),[1 3 2 4]);', obj.trafo.wavWeight));
            else
                image = permute(squeeze(bart('pics -l2', kSpaceL, ESP)),[1 3 2 4]);
            end
        else
            % create temporary dir if not yet existing
            if(~exist([obj.path,filesep,'tmpESPIRiT'],'dir'))
                mkdir([obj.path,filesep,'tmpESPIRiT']);
            end
            fftdims = size(kSpaceL);
            fftdims = fftdims(1:3);
            if(any(mod(fftdims,2) ~= 0))
                writecfl([obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], zpad(kSpaceL,size(kSpaceL) + [mod(fftdims,2),0])); % k_y - k_z - k_x - cha
            else
                writecfl([obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], kSpaceL); % k_y - k_z - k_x - cha
            end
            obj.calibSize = obj.calibSize - mod(obj.calibSize,2);
            eval(sprintf('! ecalib -r %d:%d:%d -k %d:%d:%d -m %d %s %s', obj.calibSize(1), obj.calibSize(2), obj.calibSize(3), obj.kernelSize(1), obj.kernelSize(2), obj.kernelSize(3), obj.n_maps, [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], [obj.path,filesep,'tmpESPIRiT',filesep,'espiritmaps']));
            ESP = [obj.path,filesep,'tmpESPIRiT',filesep,'espiritmaps'];

            if(strcmp(obj.solver,'L1'))
                eval(sprintf('! sense -l1 -r %f %s %s %s', obj.trafo.wavWeight, [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], ESP, [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']));
            else
                eval(sprintf('! sense -l2 %s %s %s', [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], ESP, [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']));
            end
            image = permute(squeeze(readcfl([obj.path,filesep,'tmpESPIRiT',filesep,'imgOut'])),[1 3 2 4]);
        end
%         % OR: 
%         kSpaceL = permute(ifftnshift(kSpaceL,3), [3 1 2 4]); % x - k_y - k_z - cha
%         writecfl([obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], kSpaceL); % x - k_y - k_z - cha
%         eval(sprintf('! ecalib -r %d:%d:%d -k %d:%d:%d -m %d %s %s', obj.calibSize(3), obj.calibSize(1), obj.calibSize(2), obj.kernelSize(3), obj.kernelSize(1), obj.kernelSize(2), obj.n_maps, [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], [obj.path,filesep,'tmpESPIRiT',filesep,'espiritmaps']));
%         ESP = [obj.path,filesep,'tmpESPIRiT',filesep,'espiritmaps'];
%          if(strcmp(obj.solver,'L1'))
%             eval(sprintf('! rsense -l1 -r %f %s %s %s', obj.trafo.wavWeight, [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], ESP, [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']));
%         else
%             eval(sprintf('! rsense -l2 %s %s %s', [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], ESP, [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']));
%         end   
%         image = permute(squeeze(readcfl([obj.path,filesep,'tmpESPIRiT',filesep,'imgOut'])), [2 1 3 4]);
        
        image = crop(image,[size(kSpaceL,1),size(kSpaceL,3),size(kSpaceL,2),obj.n_maps]);
        imageCha(:,1:obj.n_maps) = shiftdim(mat2cell(image,nPha,nFreq,nZ,ones(1,obj.n_maps)),2);
        [imageCha{:,obj.n_maps+1:end}] = deal(zeros(size(imageCha{1,1})));
    end
            
elseif(strcmp(obj.measPara.dimension,'4D'))
    
    % prepare k-space 
    % => k_x - k_y - k_z - cha - ... - t
    kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-3));
    kSpaceL = permute(kSpaceL,[2 1 3 5 4]);
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    % same undersampling for all coils and fully sampled in k_x direction
    obj.fullMask = abs(kSpaceL) > 0; % k_x - k_y - k_z - cha - t
    if(obj.measPara.oversampling{2,1})
        kSpaceL = ifftnshift(kSpaceL, 1); 
        kSpaceL = kSpaceL(obj.measPara.oversampling{1,1}, :, :, :, :); % anti-aliasing
        obj.fullMask = obj.fullMask(obj.measPara.oversampling{1,1}, :, :, :, :);
        kSpaceL = fftnshift(kSpaceL,1);
        kSpaceL = kSpaceL .* obj.fullMask;
        obj.measPara.dim(2) = size(kSpaceL,1);
    end
    
    if(~isempty(obj.calibSize))
        % adjust calibSize due to anti-aliasing
        obj.calibSize(3) = size(kSpaceL,1);
    end
    
    lEnsureCenter = false;
    calibSize = zeros(size(obj.fullMask,5),3);
    for iTime = 1:size(obj.fullMask,5)
        calibSize(iTime,:) = obj.calibrationSize(obj.fullMask(:,:,:,1,iTime));
        if(any(calibSize(iTime,:) <= 2))
            lEnsureCenter = true;
        end
    end
    
    if(lEnsureCenter) % ensure a center region is found!
        dPercentage = [0.065,0.2];
        m = [size(obj.fullMask,1), size(obj.fullMask,2), size(obj.fullMask,3)];
        s = [size(obj.fullMask,1)/2, round(dPercentage(1)*size(obj.fullMask,2)), round(dPercentage(2)*size(obj.fullMask,3))]; 
        idx = cell(1,length(s));
        for n=1:length(s)
            if(mod(s(n),2) == 0)
                idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
            else
                idx{n} = floor(m(n)/2)+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2)-1;
            end
        end
        lMaskNew = false(size(obj.fullMask));
        lMaskNew(idx{1},idx{2},idx{3},:,:) = true;
        lMaskNew = lMaskNew & ~obj.fullMask;
        kSpaceL(lMaskNew) = eps; 
    end
            

    kSpaceL = permute(kSpaceL,[1 2 3 4 6 7 8 9 10 11 5]);
    
    maps = bart(sprintf('ecalib -c0. -m%d',obj.n_maps), kSpaceL);
        
    bartcmd = sprintf('pics -d 5 -R W:7:0:%f',obj.lambda);
    if(obj.lambdaMC > 0)
        bartcmd = [bartcmd, sprintf(' -R T:1024:1024:%f -u 0.1', obj.lambdaMC)];
    %         bartcmd = [bartcmd, sprintf(' -R L:1024:1024:%f', obj.lambdaMC)];
    end
    if(obj.lambdaTV > 0)
        bartcmd = [bartcmd, sprintf(' -R T:7:0:%f -u 0.25', obj.lambdaTV)];
    end
    imgOut = bart(bartcmd, kSpaceL, maps);
    imgOut = permute(imgOut,[2 1 3 11 5 4 6 7 8 9 10]);
    
    imageCha{1} = imgOut;
%     imageCha = shiftdim(mat2cell(imgOut,size(imgOut,1),size(imgOut,2),size(imgOut,3),size(imgOut,4),ones(1,size(imgOut,5))),3);
end

end