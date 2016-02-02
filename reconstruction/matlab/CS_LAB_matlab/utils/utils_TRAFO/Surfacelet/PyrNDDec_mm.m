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
%%	PyrNDDec_mm.m
%%	
%%	First created: 10-11-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subs = PyrNDDec_mm(X, OutD, L, Pyr_mode, smooth_func)

%   N-dimensional multiscale pyramid decomposition - with multiple modes
%
%   INPUT:
%
%   X: input signal in the spatial domain
%  
%   OutD: "S" if we want the output to be a usual signal in the spatial "F" if
%   we are want it to be in the Fourier domain. 
%   We use this to get rid of the unnecessary (and time-consuming)
%   "ifftn-fftn" operations in the middle steps.
%
%   L: level of decomposition
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
%   subs: an L+1 by 1 cell array storing subbands from the finest scale to
%   the coarsest scale.


OutD = upper(OutD);

N = ndims(X);

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

X = fftn(X);
X = ccsym(X, N, 'c');

subs = cell(L+1, 1);

for n = 1 : L
    
    [Lp, Hp] = PrySDdec_onestep(X, w_array(n), tbw_array(n), D_array(n), smooth_func);
    
    X = Lp;
    
    Hp = ccsym(Hp, N, 'e'); 
    if OutD == 'S'
        subs{n} = real(ifftn(Hp));
    else
        subs{n} = Hp;
    end
    
    clear Lp Hp;
    
end

X = ccsym(X, N, 'e');
if OutD == 'S'
    subs{L + 1} = real(ifftn(X));
else
    subs{L + 1} = X;
end

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.       
        