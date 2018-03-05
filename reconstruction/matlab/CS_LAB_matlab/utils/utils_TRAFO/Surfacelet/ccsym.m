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
%%	ccsym.m
%%	
%%	First created: 08-14-05
%%	Last modified: 04-13-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = ccsym(x, k, type)

%% Exploit the complex conjugate symmetry in the fft of real signals.

%% type: 'c' compact 'e' expand
%% k: along which dimension

%% Dimensions of the problem
N = ndims(x);
szX = size(x);

if type == 'c'
    %% initialize the subscript array
    sub_array = repmat({':'}, [N, 1]);  
    sub_array{k} = 1 : szX(k) / 2 + 1;
    y = x(sub_array{:});
else
    %% subscript mapping for complex conjugate symmetric signal
    %% recovery
    
    szX(k) = (szX(k)-1) * 2;
    sub_conj = cell(N, 1);

    for m = 1 : N
        sub_conj{m} = [1 szX(m):-1:2];
    end
    
    sub_conj{k} = [szX(k)/2 : -1 : 2];
    %% recover the full signal in the spatial domain complex
    %% conjugate symmetry.
    y = cat(k, x, conj(x(sub_conj{:})));

end

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.
