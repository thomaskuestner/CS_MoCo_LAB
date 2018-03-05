function [imageCha, ssim_map, metrics] = algo_BFCSA_proxA_2D_real( obj,input )
% 2D real-valued variant of FCSA with an outer Bregman Iteration
% based on Huang et al. paper on FCSA
%
% input:
% obj           CS reconstruction object (holding all parameters)
% input         struct containing recon parameters and image
%
% output:
% imageCha      reconstructed channel individual image
% ssim_map      structural similarity map
% metrics       evaluation metrics 
%
% (c) Marc Fischer, Thomas Kuestner, May 2015
% -------------------------------------------------------------------------

%%
timer_proxA = tic;

%% variables:
% internal flags
flag_wavetree = true;
flag_fast = true;
flag_extendImage = true;

% internal variables:
L = 1.5; % > 1 suggested, otherwise it may diverge depending on "b{1,j} - FFT2D_mask(z{1,j},mask)"
t_old = 1;
NLTV_struct.kernelratio = 3;
NLTV_struct.windowratio = 6;
NLTV_struct.nThreads = 1;  % mind this option if used with -singleCompThread on BWCluster
itrNLTV = obj.maxitr - 5;
% chambolle tv:
parsin.MAXITER=100; parsin.tv='iso'; % 'iso' or 'l1'

% from obj:
% maxitr = obj.maxitr;
maxitr = obj.iNINNER;
% maxitrOUTER = obj.maxitrOUTER;
maxitrOUTER = obj.iNOUTER;
maxitrINNER = ceil( maxitr / maxitrOUTER );
% n1 = obj.measPara.dim(1);
% n2 = obj.measPara.dim(2);
n1 = input.n1;
n2 = input.n2;
% nSlices = obj.measPara.dim(3);
nCha = obj.measPara.dim(5);
lambdaWave = obj.lambda;
lambdaTV = obj.lambdaTV;
lambdaGroup = obj.lambdaGroup;
NLTV_struct.filterstrength = obj.lambdaNLTV_h; % 0.03 % converted from NLTV h = 0.01 %old: used: 3e-10
lambdaNLTV = obj.lambdaNLTV;
lambdaNLTV_h = obj.lambdaNLTV_h;
regularizerWeights = obj.regularizerWeights;
flagTV = obj.flagTV;
flagTV_iso = obj.flagTV_iso;
flagWave = obj.flagWave;
flagGroup = obj.flagGroup;
flagNLTV = obj.flagNLTV;
flagSBNLTV = obj.flagSBNLTV;
waveletStages = obj.trafo.waveletStages;
waveletFilterName_l1 = obj.trafo.waveletFilterName_l1;
waveletFilterName_l12 = obj.trafo.waveletFilterName_l12;

% from input:
b=input.b;
mask = input.mask;
G_prox = input.G_prox;
Gt_prox = input.Gt_prox;
groupnorm_index = input.groupnorm_index;
waveS_l1 = input.waveS_l1;
waveS_l12 = input.waveS_l12;
waveS_l12proxA = input.waveS_l12proxA;
proxA_extend_y = waveS_l12proxA(waveletStages+2,1) - waveS_l12(waveletStages+2,1);
proxA_extend_x = waveS_l12proxA(waveletStages+2,2) - waveS_l12(waveletStages+2,2);
z_proxA = cell(1,nCha);

% im_ref = input.im_ref;
% im_ref_full = zeros(n1,n2);
% for j = 1:nCha
%     im_ref_full = im_ref_full + abs(im_ref{1,j}).^2;
% end;
% im_ref_full = sqrt(im_ref_full);
clear input

% initialize cells/vectors
for j=1:nCha
    FTb{1,j} = real(iFFT2D(b{1,j}));
    % y{1,j} = real(FTb{1,j});
    % y{1,j+nCha} = imag(FTb{1,j});
end;
z = FTb; % starting point
b_outer = b;

x_wave = cell(1,nCha);
x_helper = x_wave;
x_wave_helper = x_wave;
x_tv = x_wave;
x_nltv = x_wave;
for j=1:nCha
    x_nltv{1,j} = 0;
end;
x_g = x_wave;
x_g_helper = x_wave;
x_g_proxA = x_wave;
y = x_wave;

%% MAD dependent lambdas:
flag_MAD = true;
if flag_MAD
    x_wavedec = cell(1,nCha);
    threshold = zeros(1:nCha);
    for j=1:nCha % 2*nCha
        x_wavedec{1,j} = wavedec2(z{1,j},waveletStages,waveletFilterName_l1); % atm only based on l1-daubechie
        x_wave_fine_scale = size(x_wavedec{1,j},2) - (3*waveS_l1(waveletStages+1,1)*waveS_l1(waveletStages+1,2));
        threshold(j) = mad(x_wavedec{1,j}(x_wave_fine_scale:end),1);
    end;
    clear x_wavedec
else
    threshold(1:nCha) = 1;
end;

threshold_wave(j) = lambdaWave * threshold(j) * 2/L;
threshold_TV(j) = lambdaTV * threshold(j) * 2/L;
threshold_group(j) = lambdaGroup * threshold(j) * 2/L;
threshold_NLTV(j) = lambdaNLTV; % * threshold(j) * 2/L;
threshold_NLTV_h(j) = lambdaNLTV_h; %*threshold(j); % adjust carefully or NLTV won't find a solution. lambdaNLTV_h should be < 0.01

%% initialize metrics:
itr = 0;
metrics.xtime(itr+1)= 0;
% [metrics, ssim_map{1,1}] = get_metrics_itr( im_ref, im_ref_full, z, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma );
ssim_map = [];

%% recon

dispProgress('Proximal Average', 0, maxitrOUTER*maxitrINNER);
itr = 0;
for itrOUTER = 1:maxitrOUTER  % total iter counter
    for itrINNER = 1:maxitrINNER
        
        t_new = (1+sqrt(1+4*t_old^2))/2;
        t_old = t_new;
        
        y_old = z;   % y_old = y for complex case
        
        %% landweber step
        for j = 1:nCha
            x_helper{1,j} = real(iFFT2D(FFT2D_mask(z{1,j},mask))) -FTb{1,j} ;
            z{1,j} = z{1,j} - real(x_helper{1,j})/L;
        end;
        
        %% l1-Wavelet
        if flagWave
            for j = 1:nCha % 2*nCha
                x_wave_helper{1,j} = wavedec2(z{1,j},waveletStages,waveletFilterName_l1);
                x_wave_helper{1,j} = softthresh_real(x_wave_helper{1,j},threshold_wave(j));
                x_wave{1,j} = waverec2(x_wave_helper{1,j},waveS_l1,waveletFilterName_l1);
            end;
        end;
        
        %% TV
        if flagTV
            if ~flagTV_iso
                for j = 1:nCha
                    x_tv{1,j} = MTV_2D(z{1,j},threshold_TV(j),n1,n2);
                end;
            else
                for j = 1:nCha
                    if (itr < 2)
                        [x_tv{1,j}, P]=denoise_TV_One((z{1,j}), threshold_TV(j),-inf,inf,[],parsin);
                    else
                        [x_tv{1,j}, P]=denoise_TV_One((z{1,j}), threshold_TV(j),-inf,inf,P,parsin);
                    end;
                end;
            end;
        end;
        
        %% NLTV
        if flagNLTV
            if itr >= itrNLTV
                if flagSBNLTV
                    for j = 1:nCha
                        x_nltv{1,j} = SB_NLTVfunc_slim_rescale(z{1,j},n1,n2, threshold_NLTV(j), threshold_NLTV_h(j) );
                    end;
                else
                    for j = 1:nCha
                        % if mod(itr,5) == 0 || itr == 1
                        x_nltv{1,j} = NLMF(z{1,j},NLTV_struct);
                        x_nltv{1,j} = (L.*z{1,j} + 2*threshold_NLTV(j)*x_nltv{1,j})./(L+2*threshold_NLTV(j));
                        % end;
                    end;
                end;
            end;
        end;
        
        %% l12-Wavelet
        if flagGroup
            for j = 1:nCha
                if flag_extendImage
                    z_proxA{1,j} = extend_image(z{1,j}, waveS_l12proxA, waveletStages, proxA_extend_y, proxA_extend_x);
                    x_g_helper{1,j} = wavedec2(z_proxA{1,j},waveletStages,waveletFilterName_l12);
                else
                    x_g_helper{1,j} = wavedec2(z{1,j},waveletStages,waveletFilterName_l12);
                end;
                
                if flag_wavetree
                    x_g_helper{2,j} = (G_prox*x_g_helper{1,j}')';
                else
                    x_g_helper{2,j} = zeros(1,size(x_g_helper{1,j},2));
                end;
            end;
            
            x_g_helper(1,:) = softthresh_proxA_cha(x_g_helper,threshold_group(j),nCha,waveS_l12proxA,groupnorm_index);
            
            for j = 1:nCha
                x_g_proxA{1,j} = waverec2(x_g_helper{1,j},waveS_l12proxA,waveletFilterName_l12);
                x_g{1,j} =  x_g_proxA{1,j}(1:end-proxA_extend_y,1:end-proxA_extend_x);
            end;
        end;
        
        %% add prox(.)
        for j = 1:nCha
            y{1,j} = zeros(n1,n2);
            if flagWave y{1,j} = y{1,j} + x_wave{1,j}.*regularizerWeights(1); end;
            if flagTV y{1,j} = y{1,j} + x_tv{1,j}.*regularizerWeights(2); end;
            if flagGroup y{1,j} = y{1,j} + x_g{1,j}.*regularizerWeights(3); end;
            if flagNLTV y{1,j} = y{1,j} + x_nltv{1,j}.*regularizerWeights(4); end;
            
            if itr < itrNLTV
                y{1,j} = y{1,j}/(flagTV.*regularizerWeights(2) + flagWave.*regularizerWeights(1) + flagGroup.*regularizerWeights(3));
            else
                y{1,j} = y{1,j}/(flagTV.*regularizerWeights(2) + flagNLTV.*regularizerWeights(4) + flagWave.*regularizerWeights(1) + flagGroup.*regularizerWeights(3));
            end;
            
            if flag_fast
                y{1,j}=y{1,j}+((t_old-1)/t_new).*(y{1,j}-y_old{1,j});
            end;
        end;
        
        %% metrics of current itr:
        itr = itr + 1;
%         disp(itr);
        dispProgress('Proximal Average', itr/(maxitrOUTER*maxitrINNER));
        
        metrics.xtime(itr+1)= toc(timer_proxA);
        for j = 1:nCha
            z{1,j} = y{1,j};
        end;
%         [metrics, ssim_map{1,2}] = get_metrics_itr( im_ref, im_ref_full, z, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma );
        
        if itr >= maxitr break; end;
    end;
    if itr >= maxitr break; end;
    for j=1:nCha
        b_outer{1,j} = b_outer{1,j} + b{1,j} - FFT2D_mask(z{1,j},mask);
        FTb{1,j} = real(iFFT2D(b_outer{1,j}));
    end;
end;
dispProgress('Proximal Average', 'Close');

imageCha = z;
for j = 1:nCha
    imageCha{1,j} = turn_image( imageCha{1,j} );
end;
% for j = 1:nCha+1
%     ssim_map{1,1}{1,j} = turn_image( ssim_map{1,1}{1,j} );
%     ssim_map{1,2}{1,j} = turn_image( ssim_map{1,2}{1,j} );
% end;

end

