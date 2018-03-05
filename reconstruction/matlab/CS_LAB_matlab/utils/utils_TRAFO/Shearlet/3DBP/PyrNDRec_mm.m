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
%%	PyrNDRec_mm.m
%%	
%%	First created: 10-11-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rec = PyrNDRec_mm(subs, InD, Pyr_mode, smooth_func)

%   N-dimensional multiscale pyramid reconstruction - with multiple modes
%
%   INPUT:
%
%   subs: an L+1 by 1 cell array storing subbands from the finest scale to
%   the coarsest scale.
%  
%   InD: "S" if input signal is in the spatial domain. "F" if
%   we are want it to be in the Fourier domain. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   Pyr_mode: Decomposition modes, including 1 (do not downsample the first
%   lowpass), 1.5 (downsample the first lowpass by 1.5), 2 (use a lowpass
%   with 1/3 PI cutoff frequency)
%
%   smooth_func: function handle to generate the smooth passband theta
%   function
%
%   OUTPUT:
%
%   rec: the reconstructed signal in the spatial domain


InD = upper(InD);

N = ndims(subs{1});

L = length(subs) - 1;

switch Pyr_mode
    case {1}
        w_array = [0.5 0.25 * ones(1, L - 1)];
        tbw_array = [1/6 1/12 * ones(1, L - 1)];
        D_array = [1 2 * ones(1, L - 1)];
        
    case {1.5}
        w_array = [0.5 3/8 * ones(1, L - 1)];
        tbw_array = [1/7 1/9 * ones(1, L - 1)];
        D_array = [1.5 2 * ones(1, L - 1)];
        
    case {2}
        w_array = 1 / 3 * ones(1, L);
        tbw_array = 1 / 7 * ones(1, L);
        D_array = 2 * ones(1, L);
        
    otherwise
        error('Unsupported Pyr mode.');
end

Lp = subs{end};
subs{end} = [];
if InD == 'S'
    Lp = fftn(Lp);
end
Lp = ccsym(Lp, N, 'c');

for n = L  : -1 : 1
    
    Hp = subs{n};
    subs{n} = [];
    
    if InD == 'S'
        Hp = fftn(Hp);
    end
    Hp = ccsym(Hp, N, 'c');
    Lp = PrySDrec_onestep(Lp, Hp, w_array(n), tbw_array(n), D_array(n), smooth_func);
     
end

clear Hp

Lp = ccsym(Lp, N, 'e');
rec = real(ifftn(Lp));

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
