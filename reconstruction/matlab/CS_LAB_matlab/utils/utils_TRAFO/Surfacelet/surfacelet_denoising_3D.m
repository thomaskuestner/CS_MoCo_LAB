%--------------------------------------------------------------------------
%	SurfBox-MATLAB (c)
%--------------------------------------------------------------------------
%
%	Yue M. Lu and Minh N. Do
%
%--------------------------------------------------------------------------
%
%	surfacelet_denoising_3D.m
%	
%	First created: 12-12-08
%	Last modified: 12-12-08
%
%--------------------------------------------------------------------------

function Xd = surfacelet_denoising_3D(Xn, Pyr_mode, sigma)

%   Denoising a 3-D volume using the surfacelet transform with hard
%   thresholding
%
%   Input:
%
%   Xn: a noisy 3-D volume
%
%   Pyr_mode: type of multiscale pyramid, corresponding to different levels
%   of redundancy.
%   1:       ~ 6.4
%   1.5:     ~ 4.0
%   2:       ~ 3.4
%
%   sigma: the standard deviation of the additive Gaussian noise
%
%   Output:
%
%   Xd: the denoised 3-D volume


% We use 4 levels of decomposition
L = 4;

% Different configurations for the directional decomposition
Level_64 = [-1 3 3; 3 -1 3; 3 3 -1]; % 3 * 64 directions
Level_16 = [-1 2 2; 2 -1 2; 2 2 -1]; % 3 * 64 directions
Level_4 =  [-1 1 1; 1 -1 1; 1 1 -1]; % 3 * 64 directions
Level_1 =  [-1 0 0; 0 -1 0; 0 0 -1]; % 3 directions (i.e. the hourglass decomposition)

% A parameter for the surfacelet filter
bo = 8;

switch Pyr_mode
    case 1
        Lev_array = {Level_64, Level_64, Level_16, Level_4};
    case 2
        Lev_array = {Level_64, Level_16, 1, 1};
    case 1.5
        Lev_array = {Level_64, Level_64, Level_16, Level_4};
    otherwise
        error('This Pyr_mode is invalid.');
end    

% Take the forward surfacelet transform
[Y, Recinfo] = Surfdec(Xn, Pyr_mode, Lev_array, 'ritf', 'bo', bo);
sz = size(Xn);
clear Xn;

% Load (if it exists) or build the subband scaling file used in thresholding
fname = ['SurfEP_' num2str(sz(1)) '_' num2str(sz(2)) '_' num2str(sz(3)) ...
    '_L' num2str(L) '_PyrMode' num2str(Pyr_mode) '_bo' num2str(bo) '.mat'];
if exist(fname) == 2
    load(fname);
else
    % Construct a new subband scaling file using the Monte Carlo method.
    % This can be a slow process ...
    nRepeat = 50; % number of repetitions
    Ep = SurfEP(sz, Lev_array, Pyr_mode, bo, nRepeat);
    save(fname, 'Ep');
end

% Hard-thresholding the surfacelet coefficients
for n = 1 : length(Y) - 1
    thresh = 3 * sigma + (n==1) * sigma;
    for dd = 1 : length(Y{n})
        for sd = 1 : length(Y{n}{dd})
            Y{n}{dd}{sd} = Y{n}{dd}{sd} .* (abs(Y{n}{dd}{sd}) > thresh * Ep{n}{dd}{sd});
        end
    end
end

clear Ep;

Xd = Surfrec(Y, Recinfo);
