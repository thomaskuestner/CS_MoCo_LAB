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
%%	PyrSDdec_onestep.m
%%	
%%	First created: 10-11-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Lp, Hp] = PrySDdec_onestep(X, w, tbw, D, smooth_func)

%   N-dimensional multiscale pyramid decomposition - One Step
%
%   INPUT:
%
%   X: input signal in the Fourier domain. Note: Since the signal is
%   real-vaued in the spatial domain, we remove the symmetric part of X 
%   along the last dimension to save computation time. This is important 
%   when X is a large 3-D array. See ccsym.m for details.
%  
%   w: The ideal cutoff frequency of the lowpass subband. w is specified as a
%   fraction of PI, e.g. 0.5, 1/3, ....
%
%   tbw: Width of the transition band, specified as a fraction of PI.
%
%   D: downsampling factor. For example, 1, 1.5, 2, ......
%   Note: D should divide the size of the FFT array.
%
%   smooth_func: function handle to generate the smooth passband theta
%   function
%
%   OUTPUT:
%
%   Lp: the lowpass subband in the Fourier domain with the symmetric part
%   removed along the last dimension.
%
%   Hp: the highpass subband in the Fourier domain with the symmetric part
%   removed along the last dimension.
%
%
%   Yue Lu and Minh N. Do
%   First Created: 10-11-05
%   Last Revision: 10-11-05
%

%% The dimension of the problem
N = ndims(X);
szX = size(X);

%% the original full size
szF = szX;
szF(end) = (szX(end) - 1) * 2;

%% Passband index arrays
pbd_array = cell(N, 1);

%% Lowpass filter (nonzero part)
szLf = 2 * ceil(szF / 2 * (w + tbw)) + 1;
szLf(end) = (szLf(end) + 1) / 2;
Lf = ones(szLf);

szR = ones(1, N); %% for resizing

for n = 1 : N
    %% current Fourier domain resolution
    nr = szF(n);
    
    szR(n) = szLf(n);
    
    %% passband
    pbd = [0 : ceil(nr / 2 * (w + tbw))] .';
    
    %% passband value
    pbd_value = realsqrt(smooth_func((pbd(end) - pbd) ./ (ceil(nr / 2 * (w + tbw)) ...
        - floor(nr / 2 * (w - tbw)))));
    
    %% See if we need to consider the symmetric part
    if n ~= N
        pbd = [pbd ; [nr - pbd(end) : nr - 1].'];
        pbd_value = [pbd_value ; flipud(pbd_value(2:end))];
    end
    
    pbd_array{n} = pbd + 1;
    pbd_value = reshape(pbd_value, szR);
    Lf = Lf .* repmat(pbd_value, szLf ./ szR);
    
    szR(n) = 1;
end

%% Get the lowpass subband
szLp = szF ./ D;
if any(fix(szLp) ~= szLp)
    error('The downsampling factor must be able to divide the size of the FFT matrix!');
end
szLp(end) = szLp(end) / 2 + 1;

pbd_array_sml = pbd_array;
for n = 1 : N - 1
    pbd = pbd_array{n};
    pbd((length(pbd) + 3) / 2 : end) = pbd((length(pbd) + 3) / 2 : end) ...
        + szLp(n) - szX(n);
    pbd_array_sml{n} = pbd;
end

Lp = repmat(complex(0), szLp);
Lp(pbd_array_sml{:}) = (Lf ./ D^(N/2)) .* X(pbd_array{:});


%% Get the highpass subband
Hp = X;
Hp(pbd_array{:}) = realsqrt(1 - realpow(Lf, 2)) .* X(pbd_array{:});

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
