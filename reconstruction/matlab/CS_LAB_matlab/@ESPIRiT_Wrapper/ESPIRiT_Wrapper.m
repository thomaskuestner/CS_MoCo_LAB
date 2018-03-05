classdef ESPIRiT_Wrapper < CSMaster
    % wrapper class for ESPIRiT reconstruction
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties      
        solver % reconstruction algorithm
        iNIterSplit % number of splitting iterations for CS part
        calibTyk % Tykhonov regularization parameter for calibration
        reconTyk % Tykhonov regularization in the reconstruction
        kernel % extracted GRAPPA kernel (k-space)
        kernelImg % extracted GRAPPA kernel (image space)
        eigThresh_k % eigenvalue threshold in kSpace
        eigThresh_im % eigenvalue threshold in image space
        n_maps % weighting of sensitivity maps with n eigenvalue maps
        splitWeight % splitting weight
        method % indicate which kernel should be used
        kCalib % from k-space extracted calibration data
        flagToolbox % either use ESPIRiT Linux Toolbox or MatLab implementation
        path % save path for temporary cfl files
    end
    
    methods
        function obj = ESPIRiT_Wrapper(solver, iNINNER, iNIterSplit, measPara, kernelSize, calibTyk, reconTyk, eigThresh_k, eigThresh_im, n_maps, splitWeight, trafo, flagToolbox, path, espresso, window, calibSize)
            obj = obj@CSMaster(1, iNINNER, 0.5, 0, measPara);
            obj.type = 'ESPIRiT_Wrapper';
            obj.solver = solver;
            obj.iNIterSplit = iNIterSplit;
            obj.measPara = measPara;
            obj.kernelSize = kernelSize;
            obj.calibTyk = calibTyk;
            obj.reconTyk = reconTyk;
            obj.eigThresh_k = eigThresh_k;
            obj.eigThresh_im = eigThresh_im;
            obj.n_maps = n_maps;
            obj.splitWeight = splitWeight;
            obj.trafo = trafo;
%             obj.wavWeight = wavWeight;
%             obj.waveletFilter = waveletFilter;
%             obj.waveletFilterSize = waveletFilterSize;
%             obj.waveletStages = waveletStages;
            obj.flagToolbox = flagToolbox;
            if(flagToolbox)
                currpath = fileparts(fileparts(mfilename('fullpath')));
                setenv('TOOLBOX_PATH', [currpath,filesep,'utils',filesep,'utils_BART']);
            end
            obj.espresso = espresso;
            obj.window = window;
            if(strcmp(path(end),filesep))
                obj.path = path(1:end-1);
            else
                obj.path = path;
            end
            if(nargin > 16)
                obj.calibSize = calibSize;
            end
        end
        
        function image = start(obj, dData, iLC, iTotalLines)
            image = obj.recon_main(dData, iLC, iTotalLines);
        end
    end
    
    
end

