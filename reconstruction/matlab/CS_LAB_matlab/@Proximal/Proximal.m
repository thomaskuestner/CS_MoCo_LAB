classdef Proximal < CSMaster
    % class for proximal averages
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties
%         kSpace % save kSpace in class to reduce memory overhead
%         type % reconstruction type
%         maxitr % maximum iterations = iNINNER
        itrNLTV % itr when NLTV regularizer is added
%         maxitrINNER % amount of iterations for inner loop => gets
%         calculated inside => obsolete
%         maxitrOUTER % amount of iterations for outer CS loop = iNOUTER
%         measPara % struct containing all measurement parameters (dim, dimension, ...)
%         fullMask % sampling mask
%         trafo % transformation parameter
%         lambdaWave % lagrangian multiplier for wavelet sparsity => gets lambda => obsolete
%         lambdaTV % lagrangian multiplier for TV recon
        lambdaGroup % threshold multiplier for Forest/Tree groups
        lambdaNLTV % inner NLTV parameter
        lambdaNLTV_h % weight of NLM Filter
        mue % tradeoff value for auxialary variables (in ADMM/ SB)
        regularizerWeights % relative weights of regularizers
        solver % e.g. WaTMRI, Bregmanized WaTMRI, WaTMRI_proxA, WaTMRI_SLEP, ADMM, Split-Bregman
        flagTV % TV regularizer on/off
        flagTV_iso % isotropic TV or else anisotropic
        flagWave % wavelet sparsity regularizer on/off
        flagGroup % group regularizer on/off
        flagNLTV % NLTV on/off
        flagSBNLTV % SBNLTV or NLR
        reconDIM % 2D (slice by slice) or 3D reconstruction of 3D Data
%         flagFixedChannel % reconstruct only one channel
        flagFixedSlice % determines if only a specific slice shall be reconstructed for the 2D of 3D Data case
%         flagArbitrarySubsampling % testcase subsampling for x1-Data
%         flagHoldMask % use mask from last recon -> used if algo should run with same mask for different props
%         flagMaskXY % 2D fake sampling
        flagRealAndImag % group real and imag for softthresh / groupthresh
        flagRealOnly % only reconstruct real valued image (test purpose only)
%         accel % acceleration in case of arbitrary subsampling
%         mask_variant % choose one out of 10 fixed masks
        flagSaveImages % define if recon images should be saved
        flagSaveStartImages % define if start images should be saved
%         K_1 % ssim param
%         K_2 % ssim param
%         W_size % ssim param
%         W_sigma % ssim param
    end
    
    methods
        function obj = Proximal(solver, iNINNER, itrNLTV, iNOUTER, measPara, trafo, lambda, lambdaTV, lambdaGroup, lambdaNLTV, lambdaNLTV_h, mue, regularizerWeights, flags, reconDIM, espresso)
            obj = obj@CSMaster(iNOUTER, iNINNER, 0.5, lambda, measPara);
            obj.type = 'Proximal';
            obj.solver = {solver};
%             obj.maxitr = maxitr;
            obj.itrNLTV = itrNLTV;
%             obj.maxitrINNER = maxitrINNER;
%             obj.maxitrOUTER = maxitrOUTER;
%             obj.measPara = measPara;   

            % transformation
            obj.trafo = trafo;
            obj.espresso = espresso;
            
%             obj.trafo.waveletFilter = waveletFilter;
%             obj.trafo.waveletFilterName_l1 = waveletFilterName_l1;
%             obj.trafo.waveletFilterName_l12 = waveletFilterName_l12;
%             obj.trafo.waveletStages = waveletStages;
            
            % regularizers
%             obj.lambdaWave = lambdaWave;
            obj.lambdaTV = lambdaTV;
            obj.lambdaGroup = lambdaGroup;
            obj.lambdaNLTV = lambdaNLTV;
            obj.lambdaNLTV_h = lambdaNLTV_h;
            obj.mue = mue;
            obj.regularizerWeights = regularizerWeights;
            
            % flags
            obj.flagTV = flags.flagTV;
            obj.flagTV_iso = flags.flagTV_iso;
            obj.flagWave = flags.flagWave;
            obj.flagGroup = flags.flagGroup;
            obj.flagNLTV = flags.flagNLTV;
            obj.flagSBNLTV = flags.flagSBNLTV;
            obj.reconDIM = reconDIM;
%             obj.flagFixedChannel = flags.flagFixedChannel;
            obj.flagFixedSlice = false; % no fixed slices
%             obj.flagArbitrarySubsampling = flagArbitrarySubsampling;
%             obj.flagHoldMask = flagHoldMask;
%             obj.flagMaskXY = flagMaskXY;
            obj.flagRealAndImag = flags.flagRealAndImag;
            obj.flagRealOnly = flags.flagRealOnly;
%             obj.accel = accel;
%             obj.mask_variant = mask_variant;
            obj.flagSaveImages = true; % set fixed => make it obsolete
            obj.flagSaveStartImages = false; % set fixed => make it obsolete
            
            % MSSIM parameter => obsolete
%             obj.K_1 = K_1;
%             obj.K_2 = K_2;
%             obj.W_size = W_size;
%             obj.W_sigma = W_sigma;
        end
        
        function image = recon(obj, iRep, iAvg)
            image = obj.recon_main(obj, iRep, iAvg);
        end
        
    end
    
    
end