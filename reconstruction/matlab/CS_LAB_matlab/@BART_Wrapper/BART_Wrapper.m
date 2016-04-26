classdef BART_Wrapper < CSMaster
    % wrapper class for BART reconstruction
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties              
        n_maps % weighting of sensitivity maps with n eigenvalue maps        
        kCalib % from k-space extracted calibration data
    end
    
    methods
        function obj = BART_Wrapper(lambda, lambdaCalib, lambdaTV, lambdaMC, measPara, kernelSize, n_maps, trafo, espresso, FFTwindow, currpath, calibSize)
            obj = obj@CSMaster(1, 0, 1, lambda, measPara);
            obj.type = 'BART_Wrapper';
            obj.lambdaCalib = lambdaCalib;
            obj.lambdaTV = lambdaTV;
            obj.lambdaMC = lambdaMC;
            obj.measPara = measPara;
            obj.kernelSize = kernelSize;
            obj.n_maps = n_maps;
            obj.trafo = trafo;
            setenv('TOOLBOX_PATH', [currpath,filesep,'utils',filesep,'utils_BART']);
            obj.espresso = espresso;
            obj.window = FFTwindow;
           
            if(nargin > 11)
                obj.calibSize = calibSize;
            end
        end
        
    end
    
    
end

