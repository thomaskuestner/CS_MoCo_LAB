classdef sparseMRI < CSMaster
    % wrapper class for sparseMRI/non-linear CG reconstruction
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties      
        l1Smooth % l1 smoothing strength
        lineSearchItnlim % backtracking line search step sizes alpha and beta
        lineSearchAlpha 
        lineSearchBeta % backtracking line search start point
        lineSearchT0 
        gradToll % NLCG stopping condition
        
    end
    
    methods
        function obj = sparseMRI(iNOUTER, iNINNER, measPara, lambda, lambdaTV, lambdaESPReSSo, p, trafo, espresso, window, l1Smooth, lineSearchItnlim, lineSearchAlpha, lineSearchBeta, lineSearchT0, gradToll)
            obj = obj@CSMaster(iNOUTER, iNINNER, 2*p, lambda, measPara); % use l_(2p) norm for weighting
            obj.type = 'sparseMRI';
            obj.lambdaTV = lambdaTV;
            obj.lambdaESPReSSo = lambdaESPReSSo;
%             obj.kernelSize = kernelSize;
%             obj.calibTyk = calibTyk;
%             obj.reconTyk = reconTyk;
%             obj.wavWeight = wavWeight;
            obj.trafo = trafo;
            obj.espresso = espresso;
            obj.window = window;
            obj.l1Smooth = l1Smooth;
            obj.lineSearchItnlim = lineSearchItnlim;
            obj.lineSearchAlpha = lineSearchAlpha;
            obj.lineSearchBeta = lineSearchBeta;
            obj.lineSearchT0 = lineSearchT0;
            obj.gradToll = gradToll;
        end
        
        function image = start(obj, dData, iLC, iTotalLines)
            image = obj.recon_main(dData, iLC, iTotalLines);
        end
    end
    
    
end

