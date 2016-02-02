function [ x_NLTV ] = SB_NLTVfunc_slim( x,n1,n2,lambdaNLTV_SB, lambdaNLTV_h, itr )
% solves the NLTV proximal operator
% using the SB-NLTV Code by Zhang et al.
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

% original files: see BOS_NLTV_v1: CS_NLTV, denoising_SBNLTV, update_weight

% lambdaNLTV_SB = 1; % lambda of SB prox op.
% lambdaNLTV_h = 0.0001; % opts.h0
% opts.h0=15; %% weight filter parameter, depends on noise and image standard variation. for example: for barbara [0, 255], h0=20; To be adapted for normalized image
   
% x = x;
% lambdaNLTV_SB = 0.000000001;
% lambdaNLTV_h = 0.000000001;
%% Initialization 
% Regularization parameters
% not needed - directly used lambdaNLTV_SB % opts.mu=10; %% regularization term scale
% not needed % opts.delta=1; %% delta<||A^*A||, since A here is a subsampled Fourier matrix, could be fixed to be 1.
% not needed % opts.nOuter=200; %% Outer Bregman iteration steps.
opts.nDenoising=5; %% 3st level: denoising/regularization level, in general could be fixed to be 10
% not needed % opts.type=1; %% 1: by bos,2:PBOS. 3:by Linearized Bregman If the without noise (or low), choose 1
% not needed % opts.bTol=10^-5; %%stopping criterion on residual std2(Fmask.*fft2(u)/N-Fp), if the noise standard variation is known, can set btol to be sigma
% not needed % opts.xTol=10^-5; %%stopping criterion, if the noise standard variation is known, can set btol to be sigma
opts.verbose=1; %% display message

% Weight parameters 
% not needed - directly used lambdaNLTV_h % opts.h0=15; %% weight filter parater, depends on noise and image standard variation. for example: for barbara [0, 255], h0=20; To be adapted for normalized image
opts.nWin=2; % standard: 2 %% patch size [2*nwin+1, 2*nwin+1]
opts.nBloc=7; % standard: 7 %% search window size [2*bloc+1, 2*bloc+1]
opts.nNeigh=10; % standard: 20 %% number of best neighbors (real neighbors size: 2*nNeigh+4)
opts.nWeightupdate=1; %0: no weight update, otherwise update steps 
% not needed % opts.denoising_type=1; %% NLTV denoising algorithm:1: Split bregman, 2: Projection in dual

%% Get start!
%% Main loop
  if (opts.nWeightupdate) %% update weight
      if(mod(itr,opts.nWeightupdate)==0) % only update all nWeightupdate steps

        wopts=update_weight(x,n1,n2,lambdaNLTV_h,opts.nWin,opts.nBloc,opts.nNeigh);          
      end
  end
  
  x_NLTV = denoising_SBNLTV(x,n1,n2,lambdaNLTV_SB,opts.nDenoising,wopts);

end

function wopts = update_weight(x,n1,n2,lambdaNLTV_h,nwin,nbloc,NbBestNeigh)
if (~exist('NbBestNeigh','var'))
NbBestNeigh = 10; % >=6
end    
% [n1,n2]=size(x);
NbNearNeigh = 4; % 0 or 4
display_messages = 0;
a = 2;
m=2*nwin+1;
w=2*nbloc+1;
NbNeigh=NbNearNeigh+NbBestNeigh;
G = fspecial('gaussian',[m m],a); % uses matlab filter
 VecGeneralParameters = [ display_messages; n1; n2; m; w; NbNearNeigh; NbBestNeigh; lambdaNLTV_h*lambdaNLTV_h;];
[W,Y] = compute_nl_weights_mex(single(x),single(G),single(VecGeneralParameters)); % cahnged thresholds
wopts.NbNeigh=NbNeigh;
wopts.m=m;
wopts.w=w;
wopts.W=W;
wopts.Y=Y;
end

function x_nltv = denoising_SBNLTV(x,n1,n2,mu,nDenoising,wopts)
% lambda here: variable conformity
% mu here: on fidelity term

%%
% min ( mu*|d|_1 + 1/2|u-v|_2^2 + lambda/2 | d - w(u) -b|_2^2 )
 % [n1,n2]=size(x);
 %mu=1./mu; % to be consistent with the code: mu is on the fidelity term in the denosing code
 mu = 1/mu; % reason: see above
 lambda = mu/4; % reason: see above
display_messages = 0;

% atm lambdaNLTV_SB_2 und mu vertauscht
 VecGeneralParameters = [ display_messages; n1; n2; wopts.m; wopts.w; wopts.NbNeigh; lambda; mu; nDenoising;];

d = zeros(n2,n1,wopts.NbNeigh);
b = zeros(n2,n1,wopts.NbNeigh);

x_nltv = SBNLTV_mex(single(x),single(d),single(b),single(x),...
    single(wopts.W),int32(wopts.Y),single(VecGeneralParameters));
end