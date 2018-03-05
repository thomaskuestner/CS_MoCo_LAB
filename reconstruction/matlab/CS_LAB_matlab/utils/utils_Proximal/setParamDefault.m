function [ para ] = setParamDefault( para )
% checks if all parameters are set.
% M. Fischer, April 2016

%%
    % solver: possible entries: {'WaTMRI','WaTMRI_own','WaTMRI_proxA','WaTMRI_SLEP','ADMM_proxA','SplitB'};
    if ~isfield(para,'solver') para.solver = {'FCSA_proxA'}; end

    if ~isfield(para,'cstype') para.cstype = 'CompressedS'; end
    if ~isfield(para,'reconDIM') para.reconDIM = '2D'; end
    if ~isfield(para,'waveletFilter') para.waveletFilter = '3D'; end % l1-sparsity in 2D or 3D domain
    if ~isfield(para,'waveletFilterName_l1') para.waveletFilterName_l1 = 'db2'; end
    if ~isfield(para,'waveletFilterName_l12') para.waveletFilterName_l12 = 'db1'; end
    if ~isfield(para,'waveletStages') para.waveletStages = 4; end
    if ~isfield(para,'maxitr') para.maxitr = 500; end
    if ~isfield(para,'itrNLTV') para.itrNLTV = 121; end % itr, when nltv is added. 
    if ~isfield(para,'maxitrINNER') para.maxitrINNER = 10; end
    if ~isfield(para,'maxitrOUTER') para.maxitrOUTER = 20; end
    
    % regularizer params:
    if ~isfield(para,'transform') para.transform = 'dwt'; end
    if ~isfield(para,'lambdaWave') para.lambdaWave = 0.005; end % e-10 x1, e-8 dicom
    if ~isfield(para,'lambdaTV') para.lambdaTV = 0.01; end
    if ~isfield(para,'lambdaGroup') para.lambdaGroup = 0.005; end
    if ~isfield(para,'lambdaNLTV') para.lambdaNLTV = 0.5; end % for SB sqrt(lambdaNLTV)? 
    if ~isfield(para,'lambdaNLTV_h') para.lambdaNLTV_h = 2.6; end % can fail, if chosen too small!
    if ~isfield(para,'mue') para.mue = 0.2; end % mu for ADMM / SB % tradeoff for WaTMRI
    if ~isfield(para,'regularizerWeights') para.regularizerWeights = [1 1 1 1]; end % relative weights of regularizers

    % regularizer flags:
    if ~isfield(para,'flags') para.flags = []; end;
    if ~isfield(para.flags,'flagTV') para.flags.flagTV = true; end
    if ~isfield(para.flags,'flagTV_iso') para.flags.flagTV_iso = true; end
    if ~isfield(para.flags,'flagWave') para.flags.flagWave = true; end
    if ~isfield(para.flags,'flagGroup') para.flags.flagGroup = true; end
    if ~isfield(para.flags,'flagNLTV') para.flags.flagNLTV = false; end
    if ~isfield(para.flags,'flagSBNLTV') para.flags.flagSBNLTV = true; end
    % adjust lambdaNLTV_h according to implemented code:
    if ~para.flags.flagSBNLTV para.lambdaNLTV_h = para.lambdaNLTV_h / 2550; end
    
    % sampling and save flags:   
    if ~isfield(para.flags,'flagRealAndImag') para.flags.flagRealAndImag = true; end
    if ~isfield(para.flags,'flagFixedChannel') para.flags.flagFixedChannel = false; end    
    if ~isfield(para.flags,'flagFixedSlice') para.flags.flagFixedSlice = true; end
    if ~isfield(para.flags,'flagArbitrarySubsampling') para.flags.flagArbitrarySubsampling = true; end
    if ~isfield(para.flags,'flagHoldMask') para.flags.flagHoldMask = false; end 
    if ~isfield(para.flags,'flagMaskXY') para.flags.flagMaskXY = true; end
    if ~isfield(para.flags,'flagRealOnly') para.flags.flagRealOnly = false; end
    if ~isfield(para.flags,'flagSaveImages') para.flags.flagSaveImages = true; end
    if ~isfield(para.flags,'flagSaveStartImages') para.flags.flagSaveStartImages = true; end
    if ~isfield(para.flags,'flagSaveObj') para.flags.flagSaveObj = true; end
    if ~isfield(para.flags,'flagStore') para.flags.flagStore = true; end
    if ~isfield(para.flags,'flagStoreCompact') para.flags.flagStoreCompact = false;
    if ~isfield(para,'reconDate') para.reconDate = datestr(date,'ddmm'); end
    if ~isfield(para,'accel') para.accel = 2; end
    if ~isfield(para,'mask_variant') para.mask_variant = 1; end % if flagHoldMask = true one out of 10 variants can be chosen.

    % ssim params
    if ~isfield(para,'metrics') para.metrics = []; end;
    if ~isfield(para.metrics,'K_1') para.metrics.K_1 = 0.05; end
    if ~isfield(para.metrics,'K_2') para.metrics.K_2 = 0.06; end
    if ~isfield(para.metrics,'W_size') para.metrics.W_size = 11; end
    if ~isfield(para.metrics,'W_sigma') para.metrics.W_sigma = 0.5; end
    
end

