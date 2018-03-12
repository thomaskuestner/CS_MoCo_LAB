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
%%	Surfrec.m
%%	
%%	First created: 09-20-05
%%	Last modified: 12-12-08
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rec = Surfrec(Y, Recinfo, varargin)

%   Surfacelet reconstruction
%
%   Input:
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
%   Recinfo: {Level, HGfname, bo, msize, beta, lambda, ...}  all the
%   necessary info needed for reconstruction.
%
%   IsAdjoint: (optional)
%   0 -- inverse transform     1 -- adjoint of the forward transform. 
%
%   Output:
%      
%   Rec: the reconstructed N-dimensional signal.
%
%   See also Surfdec.m
%


if (nargin == 3)
    %% Implement the adjoint operator of the forward transform. Note:
    %% NO perfect reconstruction in this case, that is, the adjoint
    %% operator is not equal to the inverse operator.
    IsAdjoint = varargin{1};
else
    %% Implement the original inverse transform.
    IsAdjoint = 0;
end


%% If possible, utilize the accelerated mex function.
%% Currently, the mex file is only compiled for Win32, Linux, and Mac.
if strcmp(Recinfo{2}{2}{4}, 'mapping')
    nomex = 1;
else
    nomex = 0;
end

if (~nomex)
    Pyr_mode = Recinfo{1};

    if (Pyr_mode == 1.5)
        Pyr_mode = 15;
    end

    Recinfo = Recinfo{2};

    L = length(Y) - 1;

    Lev_array = cell(L, 1);

    for i = 1 : L
        Lev_array{i} = int32(Recinfo{i}{3});
    end

    bo = Recinfo{1}{5};

    dec_flts = get_cbd_filters_load(bo, 'd');
    rec_flts = get_cbd_filters_load(bo, 'r');

    mSize = Recinfo{1}{6};
    beta = Recinfo{1}{7};
    lambda = Recinfo{1}{8};
    
    %% call the mex function
    Rec = mexSurfaceletRec(Y, Pyr_mode, Lev_array, dec_flts, rec_flts, mSize, beta, lambda, IsAdjoint);

    return;
end

disp(' ');
disp('Note: Pure Matlab implementation of the surfacelet transform is currently used.');
disp('Try using the accelerated mex functions instead.');
disp(' ');


Pyr_mode = Recinfo{1};
Recinfo = Recinfo{2};

L = length(Y) - 1;

%% See if we need to perform the dual-tree wavelet transform at coarse
%% scales.
L_DT = 0; %% levels of decomposition that use the dual-tree decomposition.
for n = 1 : L
    if isscalar(Recinfo{n})
        %% start using the dual-tree wavelet transform
        L_DT = n;
        break;
    end
end

if L_DT ~= 0
    Y = {Y{1:L_DT-1}  DTrec(Y(L_DT : end))} .';
    Recinfo = Recinfo(1 : L_DT - 1);
end

L = length(Y) - 1;

smooth_func = @rcos;

subs = cell(L+1, 1);

for n = L :-1: 1
    if (Recinfo{n}{1} == 'P') & (Recinfo{n}{2} == 'H')
        subs{n} = fftn(Y{n}{1}{1});
    else
        subs{n} = NDFBrec(Y{n}, 'F', Recinfo{n}, IsAdjoint);
    end
    Y{n} = [];
end

subs{L+1} = fftn(Y{L+1});

Rec = PyrNDRec_mm(subs, 'F', Pyr_mode, smooth_func);

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
