classdef SPIRiT_Wrapper < CSMaster
    % wrapper class for SPIRiT reconstruction
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties      
%         type % reconstruction type
        solver % reconstruction algorithm
%         iNINNER % amount of iterations for SPIRiT loop
%         oversampling % cut-out region in x-direction
%         dimension % 2D (x-y or x-y-t), 3D (x,y,z) or 4D (x,y,z,t) problem
%         dim % nPha x nFreq x nZ x nTime x nCha
%         measPara % struct containing all measurement parameters (dim, dimension, ...)
%         calibSize % k-space center needed for GRAPPA calibration
%         kernelSize % kernel window dimension for convolution with calibration data
        calibTyk % Tykhonov regularization parameter for calibration
        reconTyk % Tykhonov regularization in the reconstruction
        kernel % extracted GRAPPA kernel (k-space)
        kernelImg % extracted GRAPPA kernel (image space)
        method % indicate which kernel should be used
        kCalib % from k-space extracted calibration data
        A % calibration matrix
    end
    
    methods
        function obj = SPIRiT_Wrapper(solver, iNINNER, measPara, lambda, kernelSize, calibTyk, reconTyk, trafo, espresso, window, calibSize)
            obj = obj@CSMaster(1, iNINNER, 0.5, lambda, measPara);
            obj.type = 'SPIRiT_Wrapper';
            obj.solver = solver;
            obj.kernelSize = kernelSize;
            obj.calibTyk = calibTyk;
            obj.reconTyk = reconTyk;
            obj.trafo = trafo;
            obj.espresso = espresso;
            obj.window = window;
            if(nargin > 9)
                obj.calibSize = calibSize;
            end
        end
        
        function image = start(obj, dData, iLC, iTotalLines)
            image = obj.recon_main(dData, iLC, iTotalLines);
        end
    end
    
    
end

