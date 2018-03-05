function [imageCha, ssim_map, metrics] = algo_FCSA_WaTMRI_2D_real( obj,input )   
% 2D variant of WaTMRI
% based on Huang et al. paper on WaTMRI
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
timer_WaT = tic;

%% variables:
% internal flags
flag_wavetree = true;
flag_fast = true;
flag_extendImage = true;

% internal variables:
mue = 0.001;
L = 1;
t_old = 1;
NLTV_struct.kernelratio = 3;
NLTV_struct.windowratio = 6;
NLTV_struct.nThreads = 1;  % mind this option if used with -singleCompThread on BWCluster
itrNLTV = obj.iNINNER - 5;
% chambolle tv:
parsin.MAXITER=100; parsin.tv='iso'; % 'iso' or 'l1'

% from obj:
maxitr = obj.iNINNER;
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
G_wat = input.G_wat;
Gt_wat = input.Gt_wat;
% groupnorm_index = input.groupnorm_index;
waveS_l1 = input.waveS_l1;
waveS_l12 = input.waveS_l12;
waveS_WaT = input.waveS_WaT;
WaT_extend_y = waveS_WaT(waveletStages+2,1) - waveS_l12(waveletStages+2,1);
WaT_extend_x = waveS_WaT(waveletStages+2,2) - waveS_l12(waveletStages+2,2);
z_WaT = cell(1,nCha);

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
y = z;

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
x_g_WaT = x_wave;

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

for j = 1:nCha
    x_max(j) = max(max(z{1,1}))*1.0;
    x_min(j) = min(min(z{1,1}))*1.0;
end;

%% initialize metrics:
    itr = 0;
   	metrics.xtime(itr+1)= 0;    
%     [metrics, ssim_map{1,1}] = get_metrics_itr( im_ref, im_ref_full, z, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma );
    ssim_map = [];
    
%% recon
dispProgress('Proximal Average', 0, maxitr);
for itr = 1:maxitr  % total iter counter        
    
    t_new = (1+sqrt(1+4*t_old^2))/2;
    t_old = t_new;
    
    y_old = z;   % y_old = y for complex case

%% WaTMRI-Step
    if itr > 1
        if flagGroup
            for j = 1:nCha
                if flag_extendImage
                    z_WaT{1,j} = extend_image(y{1,j}, waveS_WaT, waveletStages, WaT_extend_y, WaT_extend_x);
                    x_g_helper{1,j} = wavedec2(z_WaT{1,j},waveletStages,waveletFilterName_l12);
                else
                    x_g_helper{1,j} = wavedec2(z{1,j},waveletStages,waveletFilterName_l12);    
                end;

                if flag_wavetree
                    x_g_helper{2,j} = (G_wat*x_g_helper{1,j}')'; 
                else
                    x_g_helper{2,j} = zeros(1,size(x_g_helper{1,j},2));
                end;
            end;

                x_g_helper = softthresh_group_2vec(x_g_helper,threshold_group(j),nCha); % threshold_wave(j)/threshold_group(j) acc. to author

            for j = 1:nCha
                x_g_WaT{1,j} = waverec2(x_g_helper{1,j}+(Gt_wat*x_g_helper{2,j}')',waveS_WaT,waveletFilterName_l12);
                x_g{1,j} =  x_g_WaT{1,j}(1:end-WaT_extend_y,1:end-WaT_extend_x);
            end;
        end;
    else
        for j = 1:nCha
            x_g{1,j} = zeros(n1,n2);
            z_WaT{1,j} = zeros(n1,n2);
        end;
    end;
    
%% landweber step + WaTMRI-part
    for j = 1:nCha
        x_helper{1,j} = real(iFFT2D(FFT2D_mask(z{1,j},mask))) -FTb{1,j} ;
        if itr > 1  % Phi'*G'*G*Phi*z{1,j}
            z_WaT{1,j} = extend_image(z{1,j}, waveS_WaT, waveletStages, WaT_extend_y, WaT_extend_x);
            z_WaT{1,j} = trafo_WaT(z_WaT{1,j},waveletStages,waveletFilterName_l12,G_wat,Gt_wat);
            z_WaT{1,j} = z_WaT{1,j}(1:end-WaT_extend_y,1:end-WaT_extend_x);
        end;
        z{1,j} = z{1,j} - (x_helper{1,j}  + mue*(z_WaT{1,j} - x_g{1,j}))/L; % mue used as tradeoff parameter.
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
                if (itr==1)
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
        else
            x_nltv{1,j} = 0;
       end;
    end;
    
%% add prox(.)  
    for j = 1:nCha
        y{1,j} = zeros(n1,n2);
        if flagWave y{1,j} = y{1,j} + x_wave{1,j}.*regularizerWeights(1); end;
        if flagTV y{1,j} = y{1,j} + x_tv{1,j}.*regularizerWeights(2); end;
        % if flagGroup y{1,j} = y{1,j} + x_g{1,j}.*regularizerWeights(3); end;           
        if flagNLTV y{1,j} = y{1,j} + x_nltv{1,j}.*regularizerWeights(4); end;
        
        if ~flagWave && ~flagTV && ~flagNLTV
            y{1,j} = z{1,j};
            y{1,j}(y{1,j} > x_max(j)) = x_max;
            y{1,j}(y{1,j} < x_min(j)) = x_min;
        else
            if itr < itrNLTV
                y{1,j} = y{1,j}/(flagTV.*regularizerWeights(2) + flagWave.*regularizerWeights(1)); % flagGroup.*regularizerWeights(3));
            else
                y{1,j} = y{1,j}/(flagTV.*regularizerWeights(2) + flagNLTV.*regularizerWeights(4) + flagWave.*regularizerWeights(1)); % + flagGroup.*regularizerWeights(3));
            end;
        end;
         
        flag_fast = true;
        if flag_fast
            z{1,j}=y{1,j}+((t_old-1)/t_new).*(y{1,j}-y_old{1,j});  
        end;       
    end;

%% metrics of current itr:
%     disp(itr);
    dispProgress('Proximal Average', itr/maxitr);
    
   	metrics.xtime(itr+1)= toc(timer_WaT);  
    for j = 1:nCha  
        z{1,j} = z{1,j};
    end;
%     [metrics, ssim_map{1,2}] = get_metrics_itr( im_ref, im_ref_full, z, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma );
    
end;
dispProgress('Proximal Average', 'Close');
    
    imageCha = z;
    for j = 1:nCha
        imageCha{1,j} = turn_image( imageCha{1,j} );
    end;
%     for j = 1:nCha+1
%         ssim_map{1,1}{1,j} = turn_image( ssim_map{1,1}{1,j} );
%         ssim_map{1,2}{1,j} = turn_image( ssim_map{1,2}{1,j} );
%     end;
   
end

