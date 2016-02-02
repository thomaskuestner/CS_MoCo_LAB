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
%%	NDFBrec.m
%%	
%%	First created: 04-20-05
%%	Last modified: 12-12-08
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rec = NDFBrec(Y, OutD, Recinfo, varargin)

%   NXDFB reconstruction
%
%   Y: N by 1 cell array. The k-th cell contains the subbands from
%   the k-th hourglass branch.
%
%   OutD: "S" if we want the output to be a usual signal in the spatial "F" if
%   we are want it to be in the Fourier domain. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   Recinfo = {InD, InF, Level, HGfname, bo, msize, beta, lambda ...}
%
%   InD: "S" if the input is a usual signal in the spatial domain; "F" if
%   we are given the Fourier transform of the input. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   InF: input data format. "B" (box) for a multidimensional array input; "C"
%   for a cell array input, with each cell containing one subband.
%
%   Level: an N by N matrix. The k-th row vector specifies the
%   subsequent decomposition levels for the k-th hourglass subband. The
%   individual elements of Level must be great than or equal to 0, except
%   the diagonal elements, which should be -1.
%
%   - HGfname: filter name for the hourglass filter bank. Supported types are
%     - 'ritf': rotational-invariant tight frame
%
%   - 'bo': the order of the checkerboard filters. Default: bo = 12
%   
%   - 'msize': size of the mapping kernel. This controls the quality of the
%   hourglass filters. Default = 15;
%   
%   - 'beta': beta value of the Kaiser window. This controls the quality of
%   the hourglass filters. Default = 2.5;
%
%   - 'lambda': the parameter lambda used in the hourglass filter design
%
%   IsAdjoint: (optional)
%   0 -- inverse transform     1 -- adjoint of the forward transform. 
%
%
%   Rec: the reconstructed N-dimensional signal.
%
%   See also NDFBdec.m
%


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preparation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get parameters

InD = Recinfo{1};
InF = Recinfo{2};
Level = Recinfo{3};
HGfname = Recinfo{4};
bo = Recinfo{5};
msize = Recinfo{6};
beta = Recinfo{7};
lambda = Recinfo{8};
if strcmp(HGfname,'mapping')
    norder = Recinfo{9};
    nstage = Recinfo{10};
end

if (nargin == 4)
    %% Implement the adjoint operator of the forward transform. Note:
    %% NO perfect reconstruction in this case, that is, the adjoint
    %% operator is not equal to the inverse operator.
    IsAdjoint = varargin{1};
else
    %% Implement the original inverse transform.
    IsAdjoint = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Processing Part %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dimension of the signal X
N = length(Level);
Subs = cell(N, 1);

if (~IsAdjoint)
    %% This is for the inverse transform.
    %% get the reconstruction checkerboard filters
    Gflt = get_cbd_filters_load(bo, 'r');
else
    %% This is for the adjoint operator.
    %% get the analysis checkerboard filters
    Gflt = get_cbd_filters_load(bo, 'd');
    %% then flip the filters, i.e. G_new[n] = G_old[-n].
    Gflt{2} = size(Gflt{1}) - Gflt{2} + [1 1];
    Gflt{1} = fliplr(flipud(Gflt{1}));
    Gflt{4} = size(Gflt{3}) - Gflt{4} + [1 1];
    Gflt{3} = fliplr(flipud(Gflt{3}));    
end

%% IRC reconstruction at pairs of dimensions
for k = 1 : N
    
    subband = Y{k};
    
    %% Convert the input to 'F' (Fourier) domain and 'B' (box) format
    if InD == 'F'
        if InF == 'C'
            %% just change the data format from 'Cell' to 'Box' array.
            subband = cell2box(subband, Level(k, :));
        else
            %% the easiest case: do nothing.
        end
    else
        %% First change the data format from 'Box' to 'Cell' array.
        if InF == 'B'
            subband = box2cell(subband, Level(k, :));
        end
        
        %% Take the fft and then remove the redundant part
        for n = 1 : length(subband)
            subband{n} = ccsym(fftn(subband{n}), k, 'c');
        end
        
        %% Then go back to the box format
        subband = cell2box(subband, Level(k, :));
    end
    
    %% Iteratively Resampled Checkerboard Filter Bank
    subband = IRCrecF(subband, k, Level(k,:), Gflt);
    
    Subs{k} = subband;

end

%% Last level: undecimated hourglass filter bank
if strcmp(HGfname,'ritf')
    Rec = HourGlassRec(Subs, 'F', OutD, HGfname, msize, beta, lambda);
elseif strcmp(HGfname,'mapping')
    Subs{1} = ccsym(Subs{1}, 1, 'e');
    Subs{2} = ccsym(Subs{2}, 2, 'e');
    Subs{3} = ccsym(Subs{3}, 3, 'e');
    Rec = hourglass3d_rec(Subs{1}, Subs{2}, Subs{3}, norder, nstage, 1);
    if OutD == 'S'
        Rec = real(ifftn(Rec));
    end
else
    error('Unrecognized filter type.');
end

if (IsAdjoint)
    Rec = Rec .* N;
end


%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.   
        
        