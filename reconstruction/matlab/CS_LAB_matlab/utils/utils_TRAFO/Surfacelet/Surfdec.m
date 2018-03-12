%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	SurfBox-MATLAB (c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	Yue Lu and Minh N. Do
%%
%%	Department of Electrical and Computer Engineering
%%	Coordinated Science Laboratory
%%	University of Illinois at Urbana-Champaign
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	Surfdec.m
%%	
%%	First created: 09-20-05
%%	Last modified: 12-12-08
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y, Recinfo] = Surfdec(X, Pyr_mode, Lev_array, HGfname, varargin)

%   Surfacelet decomposition
%
%   Input:
%
%   X: N-dimensional (N >= 2) input signal
%   
%   Pyr_mode: type of multiscale pyramid, corresponding to different levels
%   of redundancy.
%   1:       ~ 6.4
%   1.5:     ~ 4.0
%   2:       ~ 3.4
%
%   Lev_array: an L by 1 cell array, with each cell being an N by N matrix
%   for NDFB decomposition, or a scaler (either 1 or 2) for dual-tree
%   wavelet decomposition.
%
%   - For matrices, the k-th row vector specifies the subsequent decomposition 
%   levels for the k-th hourglass subband. The individual elements 
%   must be great than or equal to 0, except the diagonal elements, which 
%   should be -1.
%   
%   - For scalers, 1 for dual-tree real wavelet transform, and 2 for
%   dual-tree complex wavelet transform.
%
%   HGfname: filter name for the hourglass filter bank. Currently the only
%   supported type is
%   - 'ritf': rotational-invariant tight frame defined in the Fourier
%   domain.
%
%   Optional inputs: fine-tuning the NXDFB transform
%   - 'bo': the order of the checkerboard filters. Default: bo = 12
%   
%   - 'msize': size of the mapping kernel. This controls the quality of the
%   hourglass filters. Default = 15;
%   
%   - 'beta': beta value of the Kaiser window. This controls the quality of
%   the hourglass filters. Default = 2.5;
%
%   - 'lambda': the parameter lambda used in the hourglass filter design.
%   Default = 4;
%
%   Output:
%
%   Y: an L+1 by 1 cell array containing the surfacelet coefficients.
%   Y{1} to Y{L} contain bandpass coefficients, with Y{1} corresponds to
%   the finest scale.
%   Y{end} contains the lowpass subband.
%
%   For the i-th scale, if NDFB is used, Y{i} contains a 3 by 1 cell
%   array. See NDFBdec.m for details. If dualtree complex wavelet is used,
%   Y{i} contains a 2 by 1 cell array, with the first for real
%   coefficients and the second for imaginary coefficients.
%
%   Recinfo: {Level, HGfname, bo, msize, beta, lambda, ...}   
%
%   See also Surfrec.m
%


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assign default values
bo = 12;
msize = 15;
beta = 2.5;
lambda = 4;
norder = 10;
nstage = 3;

%% Read optional input arguments
k = 1; nvarargin = nargin - 4;
while k <= nvarargin
    switch varargin{k}
        case 'bo'
            if k == nvarargin
                error('Not enough input parameters.');
            else
                bo = varargin{k+1};
                k = k + 2;
            end
        case 'msize'
            if k == nvarargin
                error('Not enough input parameters.');
            else
                msize = varargin{k+1};
                k = k + 2;
            end
        case 'beta'
            if k == nvarargin
                error('Not enough input parameters.');
            else
                beta = varargin{k+1};
                k = k + 2;
            end
        case 'lambda'
            if k == nvarargin
                error('Not enough input parameters.');
            else
                lambda = varargin{k+1};
                if rem(lambda, 2) ~= 0
                    error('Lambda should be an even interger.');
                end
                k = k + 2;
            end
         case 'norder'
            if k == nvarargin
                error('Not enough input parameters.');
            else
                norder = varargin{k+1};
                k = k + 2;
            end
        case 'nstage'
            if k == nvarargin
                error('Not enough input parameters.');
            else
                nstage = varargin{k+1};
                k = k + 2;
            end
        otherwise
            error('Unrecognized input.');
    end
end

%% If possible, utilize the accelerated mex function.
%% Currently, the mex file is only compiled for Win32.

if strcmp(HGfname, 'mapping')
    nomex = 1;
else
    nomex = 0;
end

if ~nomex
    %% total level of decomposition
    L = length(Lev_array);

    %% output variables
    Recinfo = cell(L, 1);

    Lev_array_int = Lev_array;

    for i = 1 : length(Lev_array_int)
        Lev_array_int{i} = int32(Lev_array{i});
        Recinfo{i} = {'S', 'C', Lev_array{i}, 'ritf', bo, msize, beta, lambda};
    end

    Recinfo = {Pyr_mode, Recinfo};

    if (Pyr_mode == 1.5)
        Pyr_mode = 15;
    end

    dec_flts = get_cbd_filters_load(bo, 'd');
    rec_flts = get_cbd_filters_load(bo, 'r');

    %% Call the mex function
    Y = mexSurfaceletDec(X, Pyr_mode, Lev_array_int, dec_flts, rec_flts, msize, beta, lambda);

    return;
end

disp(' ');
disp('Note: Pure Matlab implementation of the surfacelet transform is currently used.');
disp('Try using the accelerated mex functions instead.');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Processing Part %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% use the raised-cosine function to get a smooth transition band. 
smooth_func = @rcos;
%% total level of decomposition
L = length(Lev_array);

%% output variables
Y = cell(L+1, 1);
Recinfo = cell(L, 1);

%% See if we need to perform the dual-tree wavelet transform at coarse
%% scales.
L_NDFB = L; %% levels of decomposition that use NDFB.
for n = 1 : L
    if isscalar(Lev_array{n})
        %% start using the dual-tree wavelet transform
        L_NDFB = n - 1;
        break;
    end
end

%% Multiscale decomposition
subs = PyrNDDec_mm(X, 'F', L_NDFB, Pyr_mode, smooth_func);
    
subs{end} = real(ifftn(subs{end}));

%% save some memory
clear X;

for n = 1 : L_NDFB
    if all(Lev_array{n} == 0)
        Y{n} = {{real(ifftn(subs{n}))}};
        subs{n} = [];
        Recinfo{n} = {'P', 'H'};
    else
        [Y_sub, Recinfo_sub] = NDFBdec(subs{n}, 'F', 'S', 'C', Lev_array{n}, HGfname, 'bo', bo, 'nstage', nstage, 'norder', norder);
        subs{n} = [];
        Y{n} = Y_sub;
        clear Y_sub;
    
        Recinfo{n} = Recinfo_sub;
    end
end

if L_NDFB == L
    %% put the lowpass in
    Y{L+1} = subs{L+1};
else
    K = L_NDFB+1;
    [Y(K : end) Recinfo(K : end)] = DTdec(subs{end}, Lev_array{K}, L - L_NDFB);
end

Recinfo = {Pyr_mode, Recinfo};

        
%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
