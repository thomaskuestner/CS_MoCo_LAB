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
%%	test_surfacelets.m
%%	
%%	First created: 04-13-06
%%	Last modified: 04-15-06
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This routine is an illustration of the basic functionality of the
%% surfacelet filter bank. Perfect reconstruction should be achieved.

N = 64;
K = 2;

%% Input data: a random 3-D array
X = randn(N, N, N);

%% Pyramid mode: 0 or 1 or 1.5 or 2, corresponding to different redundancy.
Pyr_mode = 1; %% the most redundant one (~= 6.4)

Level_H1 = [-1 K K; K -1 K; K K -1];
Level_H2 = [-1 K K; K -1 K; K K -1];

bo = 8;

%% Number of directional subbands at each scale:
%% From the finest: 3* 16, 3 * 16, ......
Lev_array = {Level_H1, Level_H1};


%% Surfacelet decomposition
tic;
[Y, Recinfo] = Surfdec(X, Pyr_mode, Lev_array, 'ritf', 'bo', bo);
toc;

%% Surfacelet reconstruction
tic
Rec = Surfrec(Y, Recinfo);
toc

%% Check perfect reconstruction
norm(X(:) - Rec(:)) / norm(X(:))

%%	This software is provided "as-is", without any express or implied
%%	warranty. In no event will the authors be held liable for any 
%%	damages arising from the use of this software.