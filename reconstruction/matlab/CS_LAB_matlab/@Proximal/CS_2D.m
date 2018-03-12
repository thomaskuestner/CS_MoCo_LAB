function [ reconData ] = CS_2D( obj, kSpace )
% prep function for 2D recon of 3D data
%
% input:
% obj           CS reconstruction object (holding all parameters)
% kSpace        subsampled kSpace
%
% output:
% reconData     cell containing reconstructed image and metrics
%
% (c) Marc Fischer, Thomas Kuestner, May 2015
% -------------------------------------------------------------------------

%% get image,b and basic parameters
nCha = obj.measPara.dim(5);

if(strcmp(obj.measPara.dimension,'3D')) % just for 3D ifft over looping direction, for 2D this would cause ifft over image dim, but in here always FFT2 applied!
    if(any(obj.trafo.fftBA(1,:)))
        for iBA=find(obj.trafo.fftBA(1,:))
            kSpace = ifftnshift(kSpace,iBA);
        end
    end
end
for j=1:nCha % choosing only a few slices doesn't work well with fft3!
    b{1,j} = kSpace(:,:,:,j); % ky - kz - kx - cha
%     mask{1,j} = obj.fullMask(:,:,:,j);
end;
[n1, n2, nSlices, tmp] = size(kSpace);

clear kSpace

%% reconstruction

iSli_start = 1;
input.mask = squeeze(obj.fullMask(:,:,1,1)); % only ky - kz mask
input.n1 = n1;
input.n2 = n2;

% struct for reconData that contains the reconCha for different solver:
nSolver = size(obj.solver,2);
solver_supported = false(1,nSolver);

for j = 1:nSolver
    %% FCSA_WaTMRI
    if(strcmp(obj.solver{1,j},'FCSA_WaTMRI'))
        % additional input variables:
        [G_wat, Gt_wat, waveS_l1, waveS_l12, waveS_WaT] = groupmatrix_ownWaT(n1,n2,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,obj.trafo.waveletFilterName_l12);
        input.G_wat = G_wat;
        input.Gt_wat = Gt_wat;
        input.waveS_l1 = waveS_l1;
        input.waveS_l12 = waveS_l12;
        input.waveS_WaT = waveS_WaT;
                
        solver_supported(j) = true;
        fprintf('Compressed Sensing reconstruction: FCSA_WaTMRI\n');
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', 0, nSlices); end;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_WaTMRI_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_WaTMRI_2D(input);
                end;
                if obj.flagSaveImages
                    FCSA_WaTMRI{1,1} = reconImageCha;
                    FCSA_WaTMRI{2,1} = reconSSIM_map;
                end;
                FCSA_WaTMRI{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_WaTMRI_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_WaTMRI_2D(input);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        FCSA_WaTMRI{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            FCSA_WaTMRI{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            FCSA_WaTMRI{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                FCSA_WaTMRI{3,nSlices +1 -s} = reconMetrics;
            end;
%             fprintf('Called the function FCSA_WaTMRI.....\n');
            if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', s/nSlices); end;
        end;
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices','Close'); end;
%         reconData.FCSA_WaTMRI = FCSA_WaTMRI;
        reconData = FCSA_WaTMRI{1,1};
        %% FCSA_SLEP
    elseif(strcmp(obj.solver{1,j},'FCSA_SLEP'))
        % additional input variables:
        [G_slep, Gt_slep, ind, nodesize, waveS_l1, waveS_l12, waveS_SLEP] = groupmatrix_SLEP(n1,n2,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,obj.trafo.waveletFilterName_l12);
        input.G_slep = G_slep;
        input.Gt_slep = Gt_slep;
        input.ind = ind;
        input.nodesize = nodesize;
        input.waveS_l1 = waveS_l1;
        input.waveS_l12 = waveS_l12;
        input.waveS_SLEP = waveS_SLEP;
        
        solver_supported(j) = true;
        fprintf('Compressed Sensing reconstruction: FCSA_SLEP\n');
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', 0, nSlices); end;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_SLEP_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_SLEP_2D(input);
                end;
                if obj.flagSaveImages
                    FCSA_SLEP{1,1} = reconImageCha;
                    FCSA_SLEP{2,1} = reconSSIM_map;
                end;
                FCSA_SLEP{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_SLEP_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_SLEP_2D(input);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        FCSA_SLEP{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            FCSA_SLEP{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            FCSA_SLEP{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                FCSA_SLEP{3,nSlices +1 -s} = reconMetrics;
            end;
%             fprintf('Called the function FCSA_SLEP.....\n');
            if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', s/nSlices); end;
        end;
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices','Close'); end;
%         reconData.FCSA_SLEP = FCSA_SLEP;
        reconData = FCSA_SLEP{1,1};
        %% FCSA_proxA
    elseif(strcmp(obj.solver{1,j},'FCSA_proxA'))
        % additional input variables:
        flag_double_nongroup_entries = true;
        [G_prox, Gt_prox, waveS_l1, waveS_l12, waveS_l12proxA, groupnorm_index] = groupmatrix_proxA(n1,n2,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,obj.trafo.waveletFilterName_l12,flag_double_nongroup_entries);
        input.groupnorm_index = groupnorm_index;
        input.G_prox = G_prox;
        input.Gt_prox = Gt_prox;
        input.waveS_l1 = waveS_l1;
        input.waveS_l12 = waveS_l12;
        input.waveS_l12proxA = waveS_l12proxA;
        
        solver_supported(j) = true;
        fprintf('Compressed Sensing reconstruction: FCSA_proxA\n');
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', 0, nSlices); end;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    FCSA_proxA{1,1} = reconImageCha;
                    FCSA_proxA{2,1} = reconSSIM_map;
                end;
                FCSA_proxA{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        FCSA_proxA{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            FCSA_proxA{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            FCSA_proxA{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                FCSA_proxA{3,nSlices +1 -s} = reconMetrics;
            end;
%             fprintf('Called the function FCSA_proxA.....\n');
            if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', s/nSlices); end;
        end;
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices','Close'); end;
%         reconData.FCSA_proxA = FCSA_proxA;
        reconData = FCSA_proxA{1,1};
        %% BFCSA_proxA
    elseif(strcmp(obj.solver{1,j},'BFCSA_proxA'))
        % additional input variables:
        flag_double_nongroup_entries = true;
        [G_prox, Gt_prox, waveS_l1, waveS_l12, waveS_l12proxA, groupnorm_index] = groupmatrix_proxA(n1,n2,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,obj.trafo.waveletFilterName_l12,flag_double_nongroup_entries);
        input.groupnorm_index = groupnorm_index;
        input.G_prox = G_prox;
        input.Gt_prox = Gt_prox;
        input.waveS_l1 = waveS_l1;
        input.waveS_l12 = waveS_l12;
        input.waveS_l12proxA = waveS_l12proxA;
        
        solver_supported(j) = true;
        fprintf('Compressed Sensing reconstruction: BFCSA_proxA\n');
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', 0, nSlices); end;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_BFCSA_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_BFCSA_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    BFCSA_proxA{1,1} = reconImageCha;
                    BFCSA_proxA{2,1} = reconSSIM_map;
                end;
                BFCSA_proxA{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_BFCSA_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_BFCSA_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        BFCSA_proxA{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            BFCSA_proxA{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            BFCSA_proxA{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                BFCSA_proxA{3,nSlices +1 -s} = reconMetrics;
            end;
%             fprintf('Called the function BFCSA_proxA.....\n');
            if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', s/nSlices); end;
        end;
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices','Close'); end;
%         reconData.BFCSA_proxA = BFCSA_proxA;
        reconData = BFCSA_proxA{1,1};
        %% ADMM_proxA
    elseif(strcmp(obj.solver{1,j},'ADMM_proxA'))
        % additional input variables:
        flag_double_nongroup_entries = true;
        [G_prox, Gt_prox, waveS_l1, waveS_l12, waveS_l12proxA, groupnorm_index] = groupmatrix_proxA(n1,n2,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,obj.trafo.waveletFilterName_l12,flag_double_nongroup_entries);
        input.groupnorm_index = groupnorm_index;
        input.G_prox = G_prox;
        input.Gt_prox = Gt_prox;
        input.waveS_l1 = waveS_l1;
        input.waveS_l12 = waveS_l12;
        input.waveS_l12proxA = waveS_l12proxA;
        
        solver_supported(j) = true;
        fprintf('Compressed Sensing reconstruction: ADMM_proxA\n');
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', 0, nSlices); end;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_ADMM_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_ADMM_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    ADMM_proxA{1,1} = reconImageCha;
                    ADMM_proxA{2,1} = reconSSIM_map;
                end;
                ADMM_proxA{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_ADMM_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_ADMM_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        ADMM_proxA{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            ADMM_proxA{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            ADMM_proxA{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                ADMM_proxA{3,nSlices +1 -s} = reconMetrics;
            end;
%             fprintf('Called the function ADMM_proxA.....\n');
            if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', s/nSlices); end;
        end;
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices','Close'); end;
%         reconData.ADMM_proxA = ADMM_proxA;
        reconData = ADMM_proxA{1,1};
        %% SB_proxA
    elseif(strcmp(obj.solver{1,j},'SB_proxA'))
        % additional input variables:
        flag_double_nongroup_entries = true;
        [G_prox, Gt_prox, waveS_l1, waveS_l12, waveS_l12proxA, groupnorm_index] = groupmatrix_proxA(n1,n2,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,obj.trafo.waveletFilterName_l12,flag_double_nongroup_entries);
        input.groupnorm_index = groupnorm_index;
        input.G_prox = G_prox;
        input.Gt_prox = Gt_prox;
        input.waveS_l1 = waveS_l1;
        input.waveS_l12 = waveS_l12;
        input.waveS_l12proxA = waveS_l12proxA;
        
        solver_supported(j) = true;
        fprintf('Compressed Sensing reconstruction: SB_proxA\n');
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', 0, nSlices); end;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
            if obj.flagFixedSlice
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_SB_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_SB_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    SB_proxA{1,1} = reconImageCha;
                    SB_proxA{2,1} = reconSSIM_map;
                end;
                SB_proxA{3,1} = reconMetrics;
            else
                for h = 1:nCha
                    input.b{1,h} = b{1,h}(:,:,s);
                end;
                if obj.flagRealOnly
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_SB_proxA_2D_real(input);
                else
                    [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_SB_proxA_2D(input);
                end;
                if obj.flagSaveImages
                    for h = 1:nCha
                        SB_proxA{1,1}{1,h}(:,:,nSlices +1 -s) = reconImageCha{1,h};
                    end;
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            SB_proxA{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            SB_proxA{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                SB_proxA{3,nSlices +1 -s} = reconMetrics;
            end;
%             fprintf('Called the function SB_proxA.....\n');
            if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices', s/nSlices); end;
        end;
        if(strcmp(obj.measPara.dimension,'3D')); dispProgress('Slices','Close'); end;
%         reconData.SB_proxA = SB_proxA;
        reconData = SB_proxA{1,1};
    %% A_PFISTA
    elseif(strcmp(obj.solver{1,j},'A_PFISTA'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l1';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_PFISTA{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_PFISTA{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_PFISTA{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_PFISTA.....\n');
        end;
%         reconData.A_PFISTA = A_PFISTA;
        reconData = A_PFISTA{1,1};
    %% A_TDIHT
    elseif(strcmp(obj.solver{1,j},'A_TDIHT'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l0';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_TDIHT{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_TDIHT{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_TDIHT{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_TDIHT.....\n');
        end;
%         reconData.A_TDIHT = A_TDIHT;
        reconData = A_TDIHT{1,1};
    %% A_ADMM_L1
    elseif(strcmp(obj.solver{1,j},'A_ADMM_L1'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l1';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_ADMM_L1{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_ADMM_L1{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_ADMM_L1{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_L1.....\n');
        end;
%         reconData.A_ADMM_L1 = A_ADMM_L1;
        reconData = A_ADMM_L1{1,1};
    %% A_ADMM_L1
    elseif(strcmp(obj.solver{1,j},'A_ADMM_L0'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'l0';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_ADMM_L0{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_ADMM_L0{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_ADMM_L0{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_L0.....\n');
        end;
%         reconData.A_ADMM_L0 = A_ADMM_L0;
        reconData = A_ADMM_L0{1,1};
    %% A_ADMM_SCAD
    elseif(strcmp(obj.solver{1,j},'A_ADMM_SCAD'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'SCAD';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        solver_supported(j) = true;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_ADMM_SCAD{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_ADMM_SCAD{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_ADMM_SCAD{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_SCAD.....\n');
        end;
%         reconData.A_ADMM_SCAD = A_ADMM_SCAD;
        reconData = A_ADMM_SCAD{1,1};
    %% A_ADMM_MCP
    elseif(strcmp(obj.solver{1,j},'A_ADMM_MCP'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'MCP';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        solver_supported(j) = true;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_ADMM_MCP{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_ADMM_MCP{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_ADMM_MCP{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_MCP.....\n');
        end;
%         reconData.A_ADMM_MCP = A_ADMM_MCP;
        reconData = A_ADMM_MCP{1,1};
    %% A_ADMM_ATAN
    elseif(strcmp(obj.solver{1,j},'A_ADMM_ATAN'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'ATAN';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        solver_supported(j) = true;
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_ADMM_ATAN{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_ADMM_ATAN{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_ADMM_ATAN{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_ATAN.....\n');
        end;
%         reconData.A_ADMM_ATAN = A_ADMM_ATAN;
        reconData = A_ADMM_ATAN{1,1};
    %% A_ADMM_PSHRINK
    elseif(strcmp(obj.solver{1,j},'A_ADMM_PSHRINK'))
        solver_supported(j) = true;
        
        % additional input variables:
        input.regularizer = 'PSHRINK';
        TransformSpecifics = get_transformSpecifics( obj.trafo.transformDict,obj.reconDIM,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,n1,n2,1);
       
        % recon
        for s = iSli_start:nSlices
%             for h = 1:nCha
%                 input.im_ref{1,h} = im_ref{1,h}(:,:,s); % only chosen slices are needed for metrics
%             end;
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
                    if(~isempty(reconSSIM_map))
                        for h = 1:nCha+1
                            A_ADMM_PSHRINK{2,1}{1,1}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,1}{1,h};
                            A_ADMM_PSHRINK{2,1}{1,2}{1,h}(:,:,nSlices +1 -s) = reconSSIM_map{1,2}{1,h};
                        end;
                    end
                end;
                A_ADMM_PSHRINK{3,nSlices +1 -s} = reconMetrics;
            end;
            fprintf('Called the function A_ADMM_PSHRINK.....\n');
        end;
%         reconData.A_ADMM_PSHRINK = A_ADMM_PSHRINK;
        reconData = A_ADMM_PSHRINK{1,1};
    end;
    if ~solver_supported(j)
        warning(['solver ' obj.solver{1,j} ' is not supported for 2D CS']);
    end;
end;

% reconData = cellfun(@(x) ipermute(x,obj.trafo.permRule), reconData, 'UniformOutput', false);


end