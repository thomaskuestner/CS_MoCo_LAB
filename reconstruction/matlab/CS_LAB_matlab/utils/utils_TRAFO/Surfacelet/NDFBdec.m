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
%%	NDFBdec.m
%%	
%%	First created: 04-18-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y, Recinfo] = NDFBdec(X, InD, OutD, OutF, Level, HGfname, varargin)

%   NDFB decomposition
%
%   X: N-dimensional (N >= 2) input signal
%
%   InD: "S" if the input X is a usual signal in the spatial domain; "F" if
%   we are given the Fourier transform of X. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   OutD: "S" if we want the output to be a usual signal in the spatial "F" if
%   we are want it to be in the Fourier domain. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   OutF: output data format. "B" for a multidimensional array output; "C"
%   for a cell array output, with each cell containing one subband.
%
%   Level: an N by N matrix. The k-th row vector specifies the
%   subsequent decomposition levels for the k-th hourglass subband. The
%   individual elements of Level must be great than or equal to 0, except
%   the diagonal elements, which should be -1.
%
%   HGfname: filter name for the hourglass filter bank. Currently supported
%   types are
%   - 'ritf': rotational-invariant tight frame  defined in the Fourier
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
%
%   Y: N by 1 cell array. The k-th (1 <= k <= N) cell is a multi-dim array,
%   containing the subbands from the k-th hourglass branch.
%
%   Recinfo: {Level, HGfname, bo, msize, beta, lambda, ...}   
%
%   See also NDFBrec.m
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


%% Read optional input arguments
k = 1; nvarargin = nargin - 6;
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
        otherwise
            error('Unrecognized input.');
    end
end

InD = upper(InD(1));
OutD = upper(OutD(1));
OutF = upper(OutF(1));

%% check the feasibility of the parameter "msize"
minsize = floor((min(size(X)) + 1) / 2);
if rem(minsize, 2) == 0
    minsize = minsize - 1;
end
if (msize > minsize)
    display(['Warning: msize changed from ' num2str(msize) ' to ' num2str(minsize) ' at the current scale.']);
    msize = minsize;
end

%% Reconstructor needs these parameters
Recinfo = {OutD, OutF, Level, HGfname, bo, msize, beta, lambda};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Processing Part %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First level: undecimated hourglass filter bank
%% Note: the output is given in the Fourier domain with the conjugate
%% symmetric part removed
Subs = HourGlassDec(X, InD, 'F', HGfname, msize, beta, lambda);

%% Dimension of the signal X
N = ndims(X);
%% Output array
Y = cell(N, 1);

   
%% get the decomposition checkerboard filters
Hflt = get_cbd_filters_load(bo, 'd');

%% Subsequent levels: IRC decomposition at pairs of dimensions
for k = 1 : N
    subband = Subs{k};

    %% free memory space
    Subs{k} = [];
    
    %% Iteratively Resampled Checkerboard Filter Bank
    subband = IRCdecF(subband, k, Level(k, :), Hflt);
        
    %% Convert the output to specified format and domain
    if OutD == 'F'
        if OutF == 'C'
            %% just change the data format from 'Box' to 'Cell' array.
            subband = box2cell(subband, Level(k, :));
        else
            %% the easiest case: do nothing.
        end
    else
        %% First change the data format from 'Box' to 'Cell' array.
        subTmp = box2cell(subband, Level(k, :));
        %% Free memory
        subband = cell(length(subTmp), 1);
        %% Perform IFFTN on each subband
        for n = 1 : length(subband)
            subband{n} = real(ifftn(ccsym(subTmp{n}, k, 'e')));
        end
        clear subTmp;
        
        if OutF == 'B'
            subband = cell2box(subband, Level(k, :));
        end
    end
    
    Y{k} = subband;

end

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
        
        