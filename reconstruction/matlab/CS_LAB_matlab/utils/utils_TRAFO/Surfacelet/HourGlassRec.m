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
%%	HourglassRec.m
%%	
%%	First created: 04-20-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rec = HourGlassRec(Subs, InD, OutD, HGfname, varargin);

%   Undecimated Hourglass Filter Bank decomposition
%
%   Subs: N by 1 cell array. The k-th cell contains either the k-th subband
%   or a string specifying the mat file name storing that subband.
%
%   InD: "S" if the input Subs are usual signals in the spatial domain; "F" if
%   we are given the Fourier transform of Subs. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   OutD: "S" if we want the output to be a usual signal in the spatial "F" if
%   we are want it to be in the Fourier domain. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   HGfname: filter name for the hourglass filter bank. Supported types are
%   - 'ritf': rotational-invariant tight frame
%
%   Optional inputs: fine-tuning the NXDFB transform
%   - 'msize': size of the mapping kernel. This controls the quality of the
%   hourglass filters. Default = 15;
%   
%   - 'beta': beta value of the Kaiser window. This controls the quality of
%   the hourglass filters. Default = 2.5;
%
%   - 'lambda': the parameter lambda used in the hourglass filter design
%
%
%   X: N-dimensional (N >= 2) input signal
%
%   See also HourGlassdec.m
%


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Praparation          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

switch HGfname
    case {'ritf', 'RITF'}
        %% Rotational-Invariant Tight Frame
        if nargin < 4
            error('Not enough input parameters');
        else
            %% Assign default values            
            msize = 15;
            beta = 2.5;
            lambda = 4;
            %% Read optional input arguments
            nvarargin = nargin - 4;
            switch nvarargin
                case {0}
                
                case {1}
                    msize = varargin{1}; 
                case {2}
                    msize = varargin{1}; 
                    beta = varargin{2}; 
                case {3}
                    msize = varargin{1}; 
                    beta = varargin{2};
                    lambda = varargin{3};
                otherwise
                    error('Unrecognized optional input argument!');
            end
        end
    otherwise
        error('Unsupported name for the hourglass filters.');
end
InD = upper(InD);
OutD = upper(OutD);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Processing Part %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dimensions of the problem
N = ndims(Subs{1});
szX = size(Subs{1});
if InD == 'F'
    %% Data in the frequency domain are half-sized.
    szX(1) = (szX(1)-1) * 2; 
end

%% For all possible 2-D diamond filters
diamond_array = cell(N-1, N-1);
%% For the N hourglass filters
hg_array = cell(N, 1);
%% denominator
denom = zeros(szX);

%% start working
switch HGfname
    case {'ritf', 'RITF'}  %% for Ratational-Invariant Tight Frame
        %% Get the Mapping Kernel
        mdiamond = NXdiamondmapping(msize, beta);
        Ldiamond = 2 * msize - 1; %% size of the diamond filter;
        %% Prepare 2-D diamond filters
        for k = 1 : N - 1
            for m = k+1 : N
                f = zeros(szX([k, m]));
                f(1:Ldiamond, 1:Ldiamond) = mdiamond;
                diamond_array{k , m - 1} = ...
                    realpow(real(fft2(circshift(f, 1 - [msize, msize]'))), lambda);
            end
        end
                
        %% We have N channels in the N-D case
        for k = 1 : N
            
            cumstrt = 0; %% a boolean switch
            
            %% work on pairs of dimensions (k,1), (k,2) .....
            for m = 1 : N
                
                if m == k
                    continue; %% except on (k,k);
                end
                
                %% Get the size of the 2-D slice
                sz = ones(1, N);
                sz([k, m]) = szX([k,m]);
                %% Get diamond filter
                p = diamond_array{min(k,m), max(k,m)-1};
                %% Make it an hourglass along dimension k
                if k < m
                    p = circshift(p, [size(p,1)/2, 0]);
                else
                    p = circshift(p, [0, size(p,2)/2]);
                end
                p = repmat(reshape(p, sz), szX ./ sz);
                if ~cumstrt
                    hg_array{k} = p;
                    cumstrt = 1;
                else
                    hg_array{k} = hg_array{k} .* p;
                end
            end
            denom = denom + hg_array{k};
        end
        
        clear p;
        
        %% we get the denominator;
        denom = N ./ denom;
               
        
        %% initialize the subscript array
        sub_array = repmat({':'}, [N, 1]);
        sub_full = sub_array;
        
        %% subscript mapping for complex conjugate symmetric signal
        %% recovery
        sub_conj = cell(N, 1);
        for k = 1 : N
                sub_conj{k} = [1 szX(k):-1:2];
        end
       
        %% get the three hourglass filters
        for k = 1 : N
            sub_array{k} = 1 : szX(k) / 2 + 1;
            hg_array{k} = ...
                realsqrt(hg_array{k}(sub_array{:}) .* denom(sub_array{:}));
            sub_array{k} = ':';
        end
        clear denom;
        
        cumstrt = 0; %% a boolean switch
        
        %% filter the signal and get the output subbands
        for k = 1 : N
            
            %% first get the subbands
            subband = Subs{k};
            Subs{k} = []; %% free the memory
            
            %% Get the nonredundant part
            sub_array{k} = 1 : szX(k) / 2 + 1;
            sub_conj{k} = [szX(k)/2 : -1 : 2];
            
            %% Convert to the Fourier domain
            if InD == 'S'
                Tmp = fftn(subband);
                subband = Tmp(sub_array{:});
                clear Tmp;
            end
               
            %% filtering ...
            Tmp = subband(sub_full{:}) .* hg_array{k};
            clear subband;
            
            hg_array{k} = [];
            
            %% grow to full size and accumulate
            if ~cumstrt
                Rec = cat(k, Tmp, conj(Tmp(sub_conj{:})));
                cumstrt = 1;
            else
                Rec = Rec + cat(k, Tmp, conj(Tmp(sub_conj{:})));
            end
            clear Tmp;
               
            sub_array{k} = ':';
            sub_conj{k} = [1 szX(k):-1:2];
        end
        
        %% If we need to get back to the spatial domain
        if OutD == 'S'
            Rec = real(ifftn(Rec)) ./ N;
        else
            Rec = Rec ./ N;
        end
            
            
    %% Currently we only have one set of filters, i.e., ritf. We are
    %% always looking for better (spatial/frequency localized) filters.
    otherwise
        error('Unsupported name for the hourglass filters.');
end

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.   