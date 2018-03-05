classdef FOCUSS < CSMaster
    % class for FOCUSS reconstruction
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties
        calibTyk % Tikhonov regularization parameter for calibration
        kernelImg % extracted GRAPPA kernel (image space)
        epsilon % noise level
        method % indicate which kernel should be used
        kCalib % from k-space extracted calibration data
        flagZeropadding % zeropad kernelImg to size(img)+kernelSize-1
        sScaling % scaling type of W matrix
        mc % Motion Correction parameters
    end
    
    methods
        function obj = FOCUSS(iNOUTER, iNINNER, p, lambda, epsilon, measPara, lambdaCalib, lambdaTV, lambdaESPReSSo, lambdaMC, kernelSize, calibTyk, calibSize, flagZeropadding, FFTwindow, trafo, espresso, sScaling, mc, lSense)
            obj = obj@CSMaster(iNOUTER, iNINNER, p, lambda, measPara);
            obj.epsilon = epsilon;
            obj.calibSize = calibSize;
            obj.lambdaCalib = lambdaCalib;
            obj.lambdaTV = lambdaTV;
            obj.lambdaESPReSSo = lambdaESPReSSo;
            obj.lambdaMC = lambdaMC;
            obj.kernelSize = kernelSize;
            obj.calibTyk = calibTyk;
            obj.flagZeropadding = flagZeropadding;
            obj.window = FFTwindow;
            obj.trafo = trafo;
            obj.espresso = espresso;
            obj.sScaling = sScaling;
            obj.mc = mc;
            obj.lSense = lSense;
            obj.type = 'FOCUSS';
        end
        
        function image = start(dData, iLC, iTotalLines, obj)
            image = recon_main(dData, iLC, iTotalLines, obj);
        end
        
        
        
    end
    
end

