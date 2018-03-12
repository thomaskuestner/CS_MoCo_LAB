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
%%	surf_coeff2vec.m
%%	
%%	First created: 04-05-07
%%	Last modified: 04-05-07
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, size_info] = surf_coeff2vec(Y)

%% Convert the coefficient of the surfacelet transform into a linear
%% vector.
%%
%% Input:
%%
%% Y: an L+1 by 1 cell array containing the surfacelet coefficients. See
%% Surfdec.m for details
%%
%%
%% Output:
%%
%% c: a linear vector storing all the coefficients.
%%
%% size_info: a book keeping cell array storing the dimension information
%% of Y. This is useful when we want to convert c back to Y. See 
%% surf_vec2coeff.m for details.
%%


%% First try to get the total number of coefficients
Ncoeff = 0;

%% number of multiscale levels
L = length(Y) - 1;
size_info = cell(L + 1, 1);

for n = 1 : L
   size_info{n} = cell(length(Y{n}), 1);
   for m = 1 : length(Y{n})
       size_info{n}{m} = cell(length(Y{n}{m}), 1); 
       for k = 1 : length(Y{n}{m})
           sz = size(Y{n}{m}{k});
           size_info{n}{m}{k} = sz;
           Ncoeff = Ncoeff + prod(sz);
       end
   end    
end

%% plus the lowpass subband
sz = size(Y{end});
size_info{L+1} = sz;
Ncoeff = Ncoeff + prod(sz);


%% Allocate the memory for the vector c
c = zeros(Ncoeff, 1);

Ncoeff = 1;

for n = 1 : L
    for m = 1 : length(Y{n})
        for k = 1 : length(Y{n}{m})
            d = prod(size(Y{n}{m}{k}));
            c(Ncoeff : Ncoeff + d - 1) = Y{n}{m}{k}(:);
            Ncoeff = Ncoeff + d;
        end
    end
end

c(Ncoeff : end) = Y{end}(:);
