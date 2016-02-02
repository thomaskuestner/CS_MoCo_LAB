function [ reconData ] = CS_3D( obj, kSpace )
% prep function for 3D recon of 3D data
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

% order kSpace in cells (one per channel)
% kSpace = permute(kSpace,obj.trafo.permRule); % sizes for sparse Dim
if(~strcmp(obj.solver{1,1},'FCSA_proxA')) % just one solver!
    if(any(obj.trafo.fftBA(1,:)))
        for iBA=find(obj.trafo.fftBA(1,:))
            kSpace = ifftnshift(kSpace,iBA);
        end
    end
end
for j=1:nCha % choosing only a few slices doesn't work well with fft3!
    b{1,j} = kSpace(:,:,:,j);
%     mask{1,j} = obj.fullMask(:,:,:,j);
end;
[n1, n2, nSlices, tmp] = size(kSpace);

clear kSpace


%% reconstruction
% struct for reconData that contains the reconCha for different solver:
nSolver = size(obj.solver,2);
solver_supported = false(1,nSolver);
iSli_start = 1;
input.n1 = n1; % sizes for sparse Dim
input.n2 = n2; % sizes for sparse Dim
input.nSlices = nSlices;
input.mask = squeeze(obj.fullMask(:,:,1,1)); % only ky - kz mask  %% MASK for sparse dim

for j = 1:nSolver
    %% FCSA_proxA
    if(strcmp(obj.solver{1,j},'FCSA_proxA'))
        % additional input variables:
        flag_double_nongroup_entries = true;
        [G_prox, Gt_prox, waveS_l1, waveS_l12, waveS_l12proxA, groupnorm_index] = groupmatrix_proxA(n1,n2,obj.trafo.waveletStages,obj.trafo.waveletFilterName_l1,obj.trafo.waveletFilterName_l12,flag_double_nongroup_entries);
        input.groupnorm_index = groupnorm_index;
        input.G_prox = G_prox;
        input.Gt_prox = Gt_prox;
        input.waveS_l1 = waveS_l1;
        input.waveS_l12 = waveS_l12;
        input.waveS_l12proxA = waveS_l12proxA;
        input.b = b;
        input.mask = squeeze(obj.fullMask(:,:,:,1)); % full 3D mask
        
        solver_supported(j) = true;
        
        fprintf('Compressed Sensing reconstruction: FCSA_proxA\n');
        [reconImageCha, reconSSIM_map, reconMetrics] = obj.algo_FCSA_proxA_3D(input);
        if obj.flagSaveImages
            FCSA_proxA{1,1} = reconImageCha;
            FCSA_proxA{2,1} = reconSSIM_map;
        end;
        FCSA_proxA{3,1} = reconMetrics;
%         fprintf('Called the function FCSA_proxA.....\n');
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
        dispProgress('Slices', 0, nSlices);
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
            dispProgress('Slices', s/nSlices);
        end;
        dispProgress('Slices','Close');
%         reconData.BFCSA_proxA = BFCSA_proxA;
        reconData = BFCSA_proxA{1,1};
%         reconData = cellfun(@(x) fftnshift(reconData,3), reconData, 'UniformOutput', false);
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
        dispProgress('Slices', 0, nSlices);
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
            dispProgress('Slices', s/nSlices);
        end;
        dispProgress('Slices','Close');
%         reconData.ADMM_proxA = ADMM_proxA;
        reconData = ADMM_proxA{1,1};
%         reconData = cellfun(@(x) fftnshift(x,3), reconData, 'UniformOutput', false);
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
        dispProgress('Slices', 0, nSlices);
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
            dispProgress('Slices', s/nSlices);
        end;
        dispProgress('Slices','Close');
%         reconData.SB_proxA = SB_proxA;
        reconData = SB_proxA{1,1};
%         reconData = cellfun(@(x) fftnshift(reconData,3), reconData, 'UniformOutput', false);
    end;
    if ~solver_supported(j)
        warning(['solver ' obj.solver{1,j} ' is not supported for 2D CS']);
    end;
end;

% reconData = cellfun(@(x) ipermute(x,obj.trafo.permRule), reconData, 'UniformOutput', false);

end