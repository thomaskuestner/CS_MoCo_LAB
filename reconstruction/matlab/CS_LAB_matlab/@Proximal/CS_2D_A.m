function [ reconData ] = CS_2D_A( obj, kSpace )
% prep function for 2D recon of 3D data

% M. Fischer, April 2016

%% get image,b and basic parameters
n1 = obj.measPara.dim(1);
n2 = obj.measPara.dim(2);
nSlices = obj.measPara.dim(3);

if ~obj.flagFixedChannel
    nCha = obj.measPara.dim(5);
else
    nCha = 1;
    obj.measPara.dim(5) = nCha;
end;

% order kSpace in cells (one per channel)
for j=1:nCha % choosing only a few slices doesn't work well with fft3!
    b{1,j} = kSpace(:,:,:,j);
    mask{1,j} = obj.fullMask(:,:,:,j);
end;

clear kSpace
obj.fullMask = [];

% choose if real image only, or complex valued one (for test purpose only)
if obj.flagRealOnly
    for j=1:nCha
        im_orig{1,j} = iFFT3D_abs(b{1,j}); % absolut value of kSpace with correct mask
    end;
else
    for j=1:nCha
        im_orig{1,j} = iFFT3D(b{1,j});
    end;
end;

%% reconstruction
% cut for faster recon:
flagCutIt = true;
if ~obj.measPara.oversampling{2,1} && flagCutIt
    if ~isempty(obj.measPara.oversampling{1,1})
        image_x_start = obj.measPara.oversampling{1,1}(1,1);
        image_x_end = obj.measPara.oversampling{1,1}(1,end);
        for j = 1:nCha
            im_orig{1,j} = im_orig{1,j}(:,image_x_start:image_x_end,:);
            b{1,j} = FFT3D(im_orig{1,j});
        end;
        n2 = image_x_end - image_x_start + 1;
        obj.measPara.dim(2) = n2;
    end;
end;

if obj.flagArbitrarySubsampling
    % arbitrary undersampling
    chosenSlice = 10;
    % chosenSlice = ceil(nSlices/2-nSlices/4); % choose an arbitrary slice of 3D mask / ceil(nSlices/2-nSlices/4)
    
    flagTestimage = false; % For testpurposes only! Shortcut to insert testimages. Overwrites chosen image
    if flagTestimage
        im_test = double(imread('MRI-Brain.jpg'));
        [n1, n2] = size(im_test);
        im_orig{1,1} = zeros(n1,n2,nSlices);
        b{1,1} = zeros(n1,n2,nSlices);
        im_orig{1,1}(:,:,chosenSlice) = im_test;
        obj.measPara.dim(1) = n1;
        obj.measPara.dim(2) = n2;
    end;
    
    
    if ~obj.flagHoldMask
        if obj.flagMaskXY
            mask = getMask_2D_2D(n1, n2, obj.accel, obj.mask_variant); % obj.mask_variant to run processes in parallel
        else
            mask = getMask_2D_1D(n1, n2, obj.accel, obj.mask_variant);
        end;
    else % throws error if file doesn't exist! use flagHoldMask = false
        if obj.flagMaskXY
            mask = getMask_2D_2D_fromTXT(n1, n2, obj.accel, obj.mask_variant); % was used with _last: now with _fromTXT to use parfor
        else
            mask = getMask_2D_1D_fromTXT(n1, n2, obj.accel, obj.mask_variant);
        end;
    end;
    
    if obj.flagFixedSlice
        for j = 1:nCha
            b{1,j}(:,:,chosenSlice) = FFT2D_mask( im_orig{1,j}(:,:,chosenSlice),mask );
            input.b{1,j} = squeeze(b{1,j}(:,:,chosenSlice));
        end;
        iSli_start = chosenSlice;
        nSlices = chosenSlice;
    else
        for j = 1:nCha
            for s = 1:nSlices
                b{1,j}(:,:,s) = FFT2D_mask( im_orig{1,j}(:,:,s),mask );
            end;
        end;
        input.b = b;
        iSli_start = 1;
    end;
else
    % has to be considered separatly -> see CS_3D. not needed atm since all test data is 3D.
    iSli_start = 1;
end;

input.mask = mask;

for j = 1:nCha % create ref and undersampled image
    im_ref{1,j} = abs(im_orig{1,j});
    for s = 1:obj.measPara.dim(3) % initial nSlices
        im_start{1,j}(:,:,s) = iFFT2D( b{1,j}(:,:,s) );
    end;
end;

if obj.flagSaveStartImages % save ref and undersampled image
    if obj.flagFixedSlice
        for j = 1:nCha
            reconData.im_start{1,j} = turn_image( im_start{1,j}(:,:,chosenSlice) );
            reconData.im_orig{1,j} = turn_image( im_orig{1,j}(:,:,chosenSlice) );
        end;
    else
        for j = 1:nCha
            reconData.im_start{1,j} = turn_image( im_start{1,j} );
            reconData.im_orig{1,j} = turn_image( im_orig{1,j} );
        end;
    end;
end;

%% recon
% struct for reconData that contains the reconCha for different solver:
nSolver = size(obj.solver,2);
solver_supported = false(1,nSolver);

for j = 1:nSolver
    %% A_PFISTA
    if(strcmp(obj.solver{1,j},'A_PFISTA'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l1';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_PFISTA{1,1} = reconImageCha;
                    A_PFISTA{2,1} = reconSSIM_map;
                end;
                A_PFISTA{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_PFISTA{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_PFISTA{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_PFISTA{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_PFISTA{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_PFISTA.....\n');
        end;
        reconData.A_PFISTA = A_PFISTA;
    %% A_TDIHT
    elseif(strcmp(obj.solver{1,j},'A_TDIHT'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l0';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_TDIHT{1,1} = reconImageCha;
                    A_TDIHT{2,1} = reconSSIM_map;
                end;
                A_TDIHT{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_PG(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_TDIHT{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_TDIHT{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_TDIHT{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_TDIHT{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_TDIHT.....\n');
        end;
        reconData.A_TDIHT = A_TDIHT;
    %% A_ADMM_L1
    elseif(strcmp(obj.solver{1,j},'A_ADMM_L1'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l1';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_ADMM_L1{1,1} = reconImageCha;
                    A_ADMM_L1{2,1} = reconSSIM_map;
                end;
                A_ADMM_L1{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_ADMM_L1{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_ADMM_L1{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_ADMM_L1{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_ADMM_L1{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_L1.....\n');
        end;
        reconData.A_ADMM_L1 = A_ADMM_L1;
    %% A_ADMM_L1
    elseif(strcmp(obj.solver{1,j},'A_ADMM_L0'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l0';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_ADMM_L0{1,1} = reconImageCha;
                    A_ADMM_L0{2,1} = reconSSIM_map;
                end;
                A_ADMM_L0{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_ADMM_L0{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_ADMM_L0{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_ADMM_L0{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_ADMM_L0{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_L0.....\n');
        end;
        reconData.A_ADMM_L0 = A_ADMM_L0;
    %% A_ADMM_SCAD
    elseif(strcmp(obj.solver{1,j},'A_ADMM_SCAD'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'SCAD';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        solver_supported(j) = true;
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_ADMM_SCAD{1,1} = reconImageCha;
                    A_ADMM_SCAD{2,1} = reconSSIM_map;
                end;
                A_ADMM_SCAD{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_ADMM_SCAD{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_ADMM_SCAD{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_ADMM_SCAD{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_ADMM_SCAD{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_SCAD.....\n');
        end;
        reconData.A_ADMM_SCAD = A_ADMM_SCAD;
    %% A_ADMM_MCP
    elseif(strcmp(obj.solver{1,j},'A_ADMM_MCP'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'MCP';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        solver_supported(j) = true;
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_ADMM_MCP{1,1} = reconImageCha;
                    A_ADMM_MCP{2,1} = reconSSIM_map;
                end;
                A_ADMM_MCP{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_ADMM_MCP{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_ADMM_MCP{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_ADMM_MCP{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_ADMM_MCP{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_MCP.....\n');
        end;
        reconData.A_ADMM_MCP = A_ADMM_MCP;
    %% A_ADMM_ATAN
    elseif(strcmp(obj.solver{1,j},'A_ADMM_ATAN'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'ATAN';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        solver_supported(j) = true;
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_ADMM_ATAN{1,1} = reconImageCha;
                    A_ADMM_ATAN{2,1} = reconSSIM_map;
                end;
                A_ADMM_ATAN{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_ADMM_ATAN{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_ADMM_ATAN{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_ADMM_ATAN{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_ADMM_ATAN{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_ATAN.....\n');
        end;
        reconData.A_ADMM_ATAN = A_ADMM_ATAN;
    %% A_ADMM_PSHRINK
    elseif(strcmp(obj.solver{1,j},'A_ADMM_PSHRINK'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'PSHRINK';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
            for h = 1:nCha
                input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
            end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    A_ADMM_PSHRINK{1,1} = reconImageCha;
                    A_ADMM_PSHRINK{2,1} = reconSSIM_map;
                end;
                A_ADMM_PSHRINK{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        A_ADMM_PSHRINK{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    for h = 1:nCha+1
                        A_ADMM_PSHRINK{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                        A_ADMM_PSHRINK{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                    end;
                end;
                A_ADMM_PSHRINK{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_PSHRINK.....\n');
        end;
        reconData.A_ADMM_PSHRINK = A_ADMM_PSHRINK;
    end;
    if ~solver_supported(j)
        fprintf(['Warning: solver ' obj.solver{1,j} ' is not supported for 2D CS']);
    end;
end;

end