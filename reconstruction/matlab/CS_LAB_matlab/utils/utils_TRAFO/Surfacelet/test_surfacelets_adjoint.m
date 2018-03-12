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
%%	test_surfacelets_adjoint.m
%%	
%%	First created: 04-06-07
%%	Last modified: 04-06-07
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This routine is an illustration of the adjoint operator of the
%% surfacelet filter bank.
%%
%% Denote by T and T* the forward surfacelet transform and its adjoint,
%% respectively. We should verify that for arbitrary vectors x and y,
%% <Tx, y>_L2 = <x, T* y>_L2.
%%

N = 64;
K1 = 2;
K2 = 1;

%% Input data: a random 3-D array
x = randn(N, N, N);

%% Pyramid mode: 0 or 1 or 1.5 or 2, corresponding to different redundancy.
Pyr_mode = 1; %% the most redundant one (~= 6.4)

Level_H1 = [-1 K1 K1; K1 -1 K1; K1 K1 -1];
Level_H2 = [-1 K2 K2; K2 -1 K2; K2 K2 -1];

bo = 8;

%% Number of directional subbands at each scale:
%% From the finest: 3* 16, 3 * 16, ......
Lev_array = {Level_H1, Level_H1};


%% Surfacelet decomposition
tic;
[x_coeff, Recinfo] = Surfdec(x, Pyr_mode, Lev_array, 'ritf', 'bo', bo);
toc;

%% Convert the coefficients to a linear vector
[x_vector, size_info] = surf_coeff2vec(x_coeff);

%% Generate the random vector y
y_vector = randn(size(x_vector));

y_coeff = surf_vec2coeff(y_vector, size_info);

%% Surfacelet reconstruction
tic
IsAdjoint = 1;
y = Surfrec(y_coeff, Recinfo, IsAdjoint);
toc

%% Check if the adjoint operator is correct
s1 = x_vector(:)' * y_vector(:);
s2 = x(:)' * y(:);

%% We are supposed to have s1 / s2 = 1.
s1 / s2


%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.