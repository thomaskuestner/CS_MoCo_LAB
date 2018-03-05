function [ reconData ] = CS_3D_A( obj, kSpace )
% prep function for 3D recon of 3D data

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
    if ~obj.flagHoldMask
        % if obj.flagMaskXY % atm no 3D_3D fct. available
        mask = getMask_3D_2D(n1, n2, obj.accel, obj.mask_variant, nSlices);
    else % throws error if file doesn't exist! use flagHoldMask = false
        % if obj.flagMaskXY
        mask = getMask_3D_2D_fromTXT(n1, n2, obj.accel, obj.mask_variant, nSlices); % does not exist atm.
    end;
    %if obj.flagFixedSlice
    for j = 1:nCha
        b{1,j} = FFT3D_mask( im_orig{1,j},mask );
    end;
    input.b = b;
else
    % nothing to do here
end;

input.mask = mask;

for j = 1:nCha % create ref and undersampled image
    im_ref{1,j} = abs(im_orig{1,j});
    im_start{1,j} = iFFT3D( b{1,j} );
end;
input.im_ref = im_ref;

if obj.flagSaveStartImages % save ref and undersampled image
    for j = 1:nCha
        reconData.im_start{1,j} = turn_image( im_start{1,j} );
        reconData.im_orig{1,j} = turn_image( im_orig{1,j} );
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
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,nSlices);
       
        % recon
        [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
        if obj.flagSaveImages
            A_PFISTA{1,1} = reconImageCha;
            A_PFISTA{2,1} = reconSSIM_map;
        end;
        A_PFISTA{3,1} = reconMetrics;
        fprintf('Called the function A_PFISTA.....\n'); 
        reconData.A_PFISTA = A_PFISTA;
    %% A_TDIHT
    elseif(strcmp(obj.solver{1,j},'A_TDIHT'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l0';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,nSlices);
       
        % recon
        [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
        if obj.flagSaveImages
            A_TDIHT{1,1} = reconImageCha;
            A_TDIHT{2,1} = reconSSIM_map;
        end;
        A_TDIHT{3,1} = reconMetrics;
        fprintf('Called the function A_TDIHT.....\n'); 
        reconData.A_TDIHT = A_TDIHT;
    %% A_ADMM_SCAD
    elseif(strcmp(obj.solver{1,j},'A_ADMM_SCAD'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'SCAD';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,nSlices);
       
        % recon
        [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
        if obj.flagSaveImages
            A_ADMM_SCAD{1,1} = reconImageCha;
            A_ADMM_SCAD{2,1} = reconSSIM_map;
        end;
        A_ADMM_SCAD{3,1} = reconMetrics;
        fprintf('Called the function A_ADMM_SCAD.....\n'); 
        reconData.A_ADMM_SCAD = A_ADMM_SCAD;
    %% A_ADMM_MCP
    elseif(strcmp(obj.solver{1,j},'A_ADMM_MCP'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'MCP';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,nSlices);
       
        % recon
        [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
        if obj.flagSaveImages
            A_ADMM_MCP{1,1} = reconImageCha;
            A_ADMM_MCP{2,1} = reconSSIM_map;
        end;
        A_ADMM_MCP{3,1} = reconMetrics;
        fprintf('Called the function A_ADMM_MCP.....\n'); 
        reconData.A_ADMM_MCP = A_ADMM_MCP;
    %% A_ADMM_ATAN
    elseif(strcmp(obj.solver{1,j},'A_ADMM_ATAN'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'ATAN';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,nSlices);
       
        % recon
        [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
        if obj.flagSaveImages
            A_ADMM_ATAN{1,1} = reconImageCha;
            A_ADMM_ATAN{2,1} = reconSSIM_map;
        end;
        A_ADMM_ATAN{3,1} = reconMetrics;
        fprintf('Called the function A_ADMM_ATAN.....\n'); 
        reconData.A_ADMM_ATAN = A_ADMM_ATAN;
    %% A_ADMM_PSHRINK
    elseif(strcmp(obj.solver{1,j},'A_ADMM_PSHRINK'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'PSHRINK';
        TransformSpecifics = get_transformSpecifics( obj.transform,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,nSlices);
       
        % recon
        [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_A_ADMM(input,TransformSpecifics);
        if obj.flagSaveImages
            A_ADMM_PSHRINK{1,1} = reconImageCha;
            A_ADMM_PSHRINK{2,1} = reconSSIM_map;
        end;
        A_ADMM_PSHRINK{3,1} = reconMetrics;
        fprintf('Called the function A_ADMM_PSHRINK.....\n'); 
        reconData.A_ADMM_PSHRINK = A_ADMM_PSHRINK;
    end;
    if ~solver_supported(j)
        fprintf(['Warning: solver ' obj.solver{1,j} ' is not supported for 2D/3D CS']);
    end;
end;

end