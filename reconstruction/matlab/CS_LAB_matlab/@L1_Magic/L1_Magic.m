classdef L1_Magic < handle
    % L1_Magic: using modified total variation of L1-magic package
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties
%         oversampling % oversampling along space directions
        type % reconstruction type
        solver % define which optimization problem should be solved
        measPara % struct containing all measurement parameters (dim, dimension, ...)
%         dimension % 2D (x-y or x-y-t), 3D (x,y,z) or 4D (x,y,z,t) problem
%         dim % nPha x nFreq x nZ x nTime x nCha  
        fullMask % mask of acquired data points
        trafo % struct of (tight frame) transformation in sparse basis
        epsilon % relaxation parameter
        lbtol % log barrier algorithm terminates when the duality gap <= lbtol
        mu % increase of barrier constant at each iteration
        cgtol % tolerance for Conjugate Gradients
        cgmaxiter % maximum number of iterations for Conjugate Gradients
        pdmaxiter % maximum number of iterations for primal dual algorithm
        pf % partial fourier parameters
        window % cut out window options
        kernelSize % kernel window dimension for convolution with calibration data
        currSlice % current slice, just needed for 2D with several slices
        currCha % current channel, just needed for 2D/3D with several channels
        lSense % apply SENSE mask for channel weighting
        kSpace % store kSpace in object
        espresso % ESPReSSo parameter
        dRecontime % recon time
    end
    
    methods
        function obj = L1_Magic(solver, measPara, epsilon, lbtol, mu, cgtol, cgmaxiter, pdmaxiter, espresso)
            obj.measPara = measPara;
            obj.epsilon = epsilon;
            obj.lbtol = lbtol;
            obj.mu = mu;
            obj.cgtol = cgtol;
            obj.cgmaxiter = cgmaxiter;  
            obj.pdmaxiter = pdmaxiter;
            obj.solver = solver;
            obj.espresso = espresso;
            obj.type = 'L1_Magic';
        end
    end
    
end

