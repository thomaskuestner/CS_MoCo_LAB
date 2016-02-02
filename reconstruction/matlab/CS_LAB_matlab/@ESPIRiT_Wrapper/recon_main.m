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
        ESP = obj.kCalib;
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
                if(~obj.flagToolbox)
                    obj.kCalib{i} = crop(kSpaceL(:,:,:,i), [obj.calibSize{i}, nCha]); % k_y - k_x - cha (cell: t)                
                end
            else
                if(~obj.flagToolbox)
                    % fixed calibration region
                    plonk = [nPha, nFreq];
                    idx = cell(1,length(plonk));
                    for iInner=1:length(plonk)
                        idx{iInner} = floor(plonk(iInner)/2)+1+ceil(-obj.calibSize(iInner)/2):floor(plonk(iInner)/2)+ceil(obj.calibSize(iInner)/2);
                    end
                    obj.kCalib{i} = kSpaceL(idx{1},idx{2},:,i);
                end
            end

            if(~obj.flagToolbox)
    %             % deprecated
    %             obj.kernel{i} = zeros([obj.kernelSize, nCha, nCha]); % k_y - k_x - cha - cha (cell: t)
    % 
    %             [AtA,obj.A{i}] = obj.corrMatrix(i);
    %             for n=1:nCha
    %                 obj.kernel{i}(:,:,:,n) = obj.calibrate(AtA,nCha,n);
    %             end

                % compute eigen-value maps
                [obj.kernel{i},S] = dat2Kernel(obj.kCalib{i},obj.kernelSize);
                idx = max(find(S >= S(1)*obj.eigThresh_k));

                % crop kernels and compute eigen-value decomposition in image space to get maps
                [M,W] = kernelEig(obj.kernel{i}(:,:,:,1:idx),[nPha,nFreq]);

                % Compute Soft-SENSE ESPIRiT Maps 
                % crop sensitivity maps according to eigenvalues == 1
                % weigth with n maps of eigen-values

                maps = M(:,:,:,end-(obj.n_maps-1):end);

                % Weight the eigenvectors with soft-senses eigen-values
                weights = W(:,:,end-(obj.n_maps-1):end) ;
                weights = (weights - obj.eigThresh_im)./(1-obj.eigThresh_im).* (W(:,:,end-(obj.n_maps-1):end) > obj.eigThresh_im);
                weights = -cos(pi*weights)/2 + 1/2;

                % create and ESPIRiT operator
                ESP{i} = ESPIRiT(maps,weights);
            else
                if(ispc)
                    eval(sprintf('ESP{i} = bart(''ecalib -r %d:%d -k %d:%d -m %d -1'',kSpaceL(:,:,:,i));', obj.calibSize{i}(1), obj.calibSize{i}(2), obj.kernelSize(1), obj.kernelSize(2), obj.n_maps ));

                else
                    % create temporary dir if not yet existing
                    if(~exist([obj.path,filesep,'tmpESPIRiT'],'dir'))
                        mkdir([obj.path,filesep,'tmpESPIRiT']);
                    end
                    writecfl([obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], kSpaceL(:,:,:,i)); % k_y - k_x - cha
                    eval(sprintf('! ecalib -r %d:%d -k %d:%d -m %d %s %s', obj.calibSize{i}(1), obj.calibSize{i}(2), obj.kernelSize(1), obj.kernelSize(2), obj.n_maps, [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], [obj.path,filesep,'tmpESPIRiT',filesep,'espiritmaps_',num2str(i)]));
                    ESP{i} = [obj.path,filesep,'tmpESPIRiT',filesep,'espiritmaps_',num2str(i)];
                end
            end
        end

        % ensure that empty entries are zero
%         mask = false(nPha,nFreq,nCha,nTime);
        mask = obj.fullMask;
%         for i=1:nTime
%             mask(:,:,:,i) = repmat(obj.fullMask(:,:,i),[1,1,nCha]); 
%         end
        kSpaceL(~mask) = 0;
                  
        if(strcmp(obj.solver,'L1') && ~obj.flagToolbox)
            if(strcmp(obj.trafo.trafoType,'wavelet_lab'))
                ssx = 2^ceil(log2(nPha)); 
                ssy = 2^ceil(log2(nFreq));
                ss = max(ssx, ssy);
                XOP = Wavelet(obj.trafo.waveletFilter,obj.trafo.waveletFilterSize,log2(ss/2^(obj.trafo.waveletStages)));
            end
        end
        
        image = zeros(nPha,nFreq,obj.n_maps,nTime);
        dispProgress('Time', 0, nTime);
        fprintf('ESPIRiT reconstruction\n');
        for i=1:nTime
            if(~obj.flagToolbox)
                if(strcmp(obj.solver,'L1'))
                    FT = p2DFT(mask(:,:,:,i),[nPha,nFreq,nCha]);
                    if(strcmp(obj.trafo.trafoType,'fft'))
                        XOP = FT;
                    end
                    image(:,:,:,i) = obj.cgL1ESPIRiT(kSpaceL(:,:,:,i), zeros(size(kSpaceL(:,:,:,i))), FT, ESP, obj.iNINNER, XOP, obj.trafo.wavWeight, obj.splitWeight, obj.iNIterSplit); % => y-x-cha-t
                else
                    [~,image(:,:,:,i)] = obj.cgESPIRiT(kSpaceL(:,:,:,i),ESP, obj.iNINNER, obj.calibTyk, zeros(size(kSpaceL(:,:,:,i))));
                end
            else
                if(ispc)
                    if(strcmp(obj.solver,'L1'))
                        eval(sprintf('image(:,:,:,i) = squeeze(bart(''sense -l1 -r %f'', kSpaceL(:,:,:,i), ESP{i}));', obj.trafo.wavWeight));
                    else
                        image(:,:,:,i) = squeeze(bart('sense -l2', kSpaceL(:,:,:,i), ESP{i}));
                    end
                else
                    fftdims = size(kSpaceL);
                    fftdims = fftdims(1:2);
                    if(any(mod(fftdims,2) ~= 0))
                        writecfl([obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], zpad(kSpaceL(:,:,:,i),size(kSpaceL(:,:,:,i)) + [mod(fftdims,2),0])); % k_y - k_z - k_x - cha
                    else
                        writecfl([obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], kSpaceL(:,:,:,i)); % k_y - k_x - cha
                    end
                    if(strcmp(obj.solver,'L1'))
                        eval(sprintf('! sense -l1 -r %f %s %s %s', obj.trafo.wavWeight, [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], ESP{i}, [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']));
                    else
                        eval(sprintf('! sense -l2 %s %s %s', [obj.path,filesep,'tmpESPIRiT',filesep,'kspace'], ESP{i}, [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']));
                    end
    %                 if(obj.n_maps > 1)
    %                     % combine ESPIRiT images into one
    %                    eval(sprintf('! rss 3 %s %s', [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut'], [obj.path,filesep,'tmpESPIRiT',filesep,'imgOutRSS']));
    %                    delete([obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']);
    %                    movefile([obj.path,filesep,'tmpESPIRiT',filesep,'imgOutRSS'], [obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']);
    %                 end
                    image(:,:,:,i) = squeeze(readcfl([obj.path,filesep,'tmpESPIRiT',filesep,'imgOut']));
                end
            end
            dispProgress('Time',i/nTime);
        end
        dispProgress('Time','Close');
        
        imageCha(iSli,1:obj.n_maps) = shiftdim(mat2cell(permute(image,[1 2 4 3]),nPha,nFreq,nTime,ones(1,obj.n_maps)),2);
        [imageCha{iSli,obj.n_maps+1:end}] = deal(zeros(size(imageCha{iSli,1}))); 

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
    % TODO

end

if(obj.flagToolbox)
    rmdir([obj.path,filesep,'tmpESPIRiT'], 's');
end

end