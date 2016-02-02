%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%	SurfBox-MATLAB (c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	Yue M. Lu and Minh N. Do
%%
%%	Department of Electrical and Computer Engineering
%%	Coordinated Science Laboratory
%%	University of Illinois at Urbana-Champaign
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%	surf_vec2coeff.m
%%	
%%	First created: 04-06-07
%%	Last modified: 04-06-07
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = surf_coeff2vec(c, size_info)

%% Convert the surfacelet coefficient stored in a linear vector to a nested
%% cell array.
%%
%% Input:
%%
%% c: a linear vector storing all the coefficients.
%%
%% size_info: a book keeping cell array storing the dimension information
%% of Y. This is useful when we want to convert c back to Y. See 
%% surf_vec2coeff.m for details.
%%
%% Output:
%%
%% Y: an L+1 by 1 cell array containing the surfacelet coefficients. See
%% Surfdec.m for details
%%


%% we can inherit the nested structure
Y = size_info;

%% number of multiscale levels
L = length(Y) - 1;

Ncoeff = 1;

for n = 1 : L
    for m = 1 : length(Y{n})
        for k = 1 : length(Y{n}{m})
            sz = size_info{n}{m}{k};
            d = prod(sz);
            Y{n}{m}{k} = reshape(c(Ncoeff : Ncoeff + d - 1), sz);
            Ncoeff = Ncoeff + d;
        end
    end
end

%% the lowpass subband
sz = size_info{end};
Y{end} = reshape(c(Ncoeff : end), sz);




