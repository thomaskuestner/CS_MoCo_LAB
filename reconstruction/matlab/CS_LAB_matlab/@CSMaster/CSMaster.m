classdef CSMaster < handle
    %CS_LAB superclass
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties
        type % reconstruction type
        kSpace % save kSpace in class to reduce memory overhead
        iNOUTER % outer loop iterations
        iNINNER % inner loop iterations
        p % L-p norm for reweighted minimum norm
        lambda % lagrange multiplier for noise stability
        lambdaCalib % lagrange multiplier for calibration
        lambdaTV % lagrange multiplier for TV
        lambdaESPReSSo % lagrange multiplier for ESPReSSo constraint
        lambdaMC % lagrange multiplier for motion correction
%         oversampling % oversampling along space directions
%         dimension % 2D (x-y or x-y-t), 3D (x,y,z) or 4D (x,y,z,t) problem
%         dim % nPha x nFreq x nZ x nTime x nCha
        measPara % struct containing all measurement parameters (dim, dimension, ...)
        kernelSize % kernel window dimension for convolution with calibration data
        calibSize % k-space center needed for initial estimate of W or GRAPPA calibration
        fullMask % mask of acquired data points
        currSlice % current slice, just needed for 2D with several slices
        currCha % current channel, just needed for 2D/3D with several channels
        trafo % struct of (tight frame) transformation in sparse basis
        espresso % espresso parameters
        window % cut out window options
        lSense % apply sense mask for multicoil reconstruction
        dRecontime % reconstruction time
    end
    
    methods
        function obj = CSMaster(iNOUTER, iNINNER, p, lambda, measPara)
            % constructor
            obj.iNOUTER = iNOUTER;
            obj.iNINNER = iNINNER;
            obj.p       = p;
            obj.lambda  = lambda;
            obj.measPara = measPara;
        end
    end
        
end

