%--------------------------------------------------------------------------
%	SurfBox-MATLAB (c)
%--------------------------------------------------------------------------
%
%	Yue M. Lu and Minh N. Do
%
%--------------------------------------------------------------------------
%
%	SurfEP.m
%	
%	First created: 12-12-08
%	Last modified: 12-12-08
%
%--------------------------------------------------------------------------

function Ep = SurfEP(sz, Lev_array, Pyr_mode, bo, nRepeat)

%   Estimating the subband scaling coefficients of the surfacelet transform
%   using a Monte Carlo method.
%
%   Input:
%
%   sz: the size of a N-dimensional array (N >= 2)
%
%   Lev_array: an L by 1 cell array, with each cell being an N by N matrix
%   for NDFB decomposition, or a scaler (either 1 or 2) for dual-tree
%   wavelet decomposition.
%
%   Pyr_mode: type of multiscale pyramid, corresponding to different levels
%   of redundancy.
%   1:       ~ 6.4
%   1.5:     ~ 4.0
%   2:       ~ 3.4
%
%   bo: the order of the checkerboard filters. Default: bo = 12
%
%   nRepeat: the number of trials in the Monte Carlo method.
%
%   Output:
%
%   Ep: the estimated subband scaling coefficients


h = waitbar(0, ['Estimating the subband scaling coefficients.' ...
        'This can be a slow process, but it only needs to be done once for a given configuration.']);

for n = 1 : nRepeat
    X = randn(sz);
    [Y, Recinfo] = Surfdec(X, Pyr_mode, Lev_array, 'ritf', 'bo', bo);
    clear X; % to save some memory
    
    if n == 1
        Ep = Y;
        
        for j = 1 : length(Y) - 1
            for dd = 1 : length(Y{j})
                for sd = 1 : length(Y{j}{dd})
                    Ep{j}{dd}{sd} = Ep{j}{dd}{sd} .^ 2;
                end
            end
        end
        
    else
        for j = 1 : length(Y) - 1
            for dd = 1 : length(Y{j})
                for sd = 1 : length(Y{j}{dd})
                    Ep{j}{dd}{sd} = Ep{j}{dd}{sd} + Y{j}{dd}{sd} .^ 2;
                end
            end
        end
    end
    
    clear Y;
    
    waitbar(n / nRepeat, h);
end
close(h);

for j = 1 : length(Ep) - 1
    for dd = 1 : length(Ep{j})
        for sd = 1 : length(Ep{j}{dd})
            Ep{j}{dd}{sd} = sqrt(Ep{j}{dd}{sd} / (nRepeat - 1));
        end
    end
end

