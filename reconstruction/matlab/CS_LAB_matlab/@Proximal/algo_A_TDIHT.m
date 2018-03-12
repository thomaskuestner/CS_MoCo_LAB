function [imageCha, ssim_map, metrics] = algo_A_TDIHT( obj,input )
%%
% 3D variant of FCSA
% based on Huang et al. paper on FCSA
% Marc Fischer, May 2015

%%
timer_proxA = tic;



%% variables:
flag_fast = false;
L = 5;
regularizer = 'l0';

% % from obj:
% maxitr = obj.maxitr;
% n1 = obj.measPara.dim(1);
% n2 = obj.measPara.dim(2);
% %% temporary CHANGE!
nSlices = 1; % obj.measPara.dim(3);
% nCha = obj.measPara.dim(5);
% lambdaWave = obj.lambdaWave;
maxitr = obj.iNINNER;
% n1 = obj.measPara.dim(1);
% n2 = obj.measPara.dim(2);
n1 = input.n1;
n2 = input.n2;
if strcmp(obj.reconDIM,'3D')
    n3 = obj.measPara.dim(2);
elseif strcmp(obj.reconDIM,'2D')
    n3 = 1;
end;
% nSlices = obj.measPara.dim(3);
nCha = obj.measPara.dim(5);
lambdaWave = obj.lambda;
tau = obj.mue;

flagRealAndImag = obj.flagRealAndImag;
waveletStages = obj.trafo.waveletStages;
maxWavecells = waveletStages*7+1;
waveletFilterName_l1 = obj.trafo.waveletFilterName_l1;
waveS_l1 = input.waveS_l1;

%% operators
ST = @(x) operator_frame(x,n1,n2,1,1,0, @(x) sparse_trans(x,waveletStages,waveletFilterName_l1));
BT = @(x) operator_frame(x,n1,n2,1,0,1, @(x) back_trans(x,waveS_l1,waveletFilterName_l1)); % CARE n1, n2 IS DEPENDENT ON COEFF NUMBERS!

%% from input:
b=input.b;
mask = input.mask;



im_ref = input.im_ref;
im_ref_full = zeros(n1,n2,nSlices);
for j = 1:nCha
    im_ref_full = im_ref_full + abs(im_ref{1,j}).^2;
end;
im_ref_full = sqrt(im_ref_full);
clear input

% initialize cells/vectors
for j=1:nCha
    At_v{1,j} = iFFT2D(b{1,j}); % !!!!!!!!!!!!!!!!!!! Temporary: adjust for correct dim!
end;

x = At_v;
x_accel = x;
x_old = x;
z_coeff = cell(1,nCha);  
AtA_x = cell(1,nCha);
x_test = cell(1,nCha);

t_new = 1;

%% MAD dependent lambdas:
flag_MAD = true;
if flag_MAD
    x_wavedec_2D = cell(1,nCha);
    x_wavedec_3D = cell(1,nCha);
    threshold_2D = zeros(nCha,1); % one threshold value is enough (only for testing purposes)
    threshold_3D = zeros(nCha,1); % one threshold value is enough (only for testing purposes)
    if flagRealAndImag
        for j=1:nCha
            for l = 1:nSlices
                x_wavedec_2D{1,j} = wavedec2(abs(At_v{1,j}),waveletStages,waveletFilterName_l1); % atm only based on l1-daubechie
                x_wave_fine_scale_2D = size(x_wavedec_2D{1,j},2) - (3*waveS_l1(waveletStages+1,1)*waveS_l1(waveletStages+1,2));
                threshold_2D(j) = mad(x_wavedec_2D{1,j}(x_wave_fine_scale_2D:end),1);
            end;
%             x_wavedec_3D{1,j} = wavedec3(abs(z{1,j}),waveletStages,waveletFilterName_l1); % atm only based on l1-daubechie
%             x_wave_fine_scale_3D = zeros(0,0);
%             for k = (7*(x_wavedec_3D{1,j}.level-1)+2):maxWavecells
%                 x_wave_fine_scale_3D = [x_wave_fine_scale_3D  x_wavedec_3D{1,j}.dec{k}];
%             end;
%             [size1, size2, size3] = size(x_wave_fine_scale_3D);
%             threshold_3D(j) = mad(reshape(x_wave_fine_scale_3D,size1*size2*size3,1));
        end;
    else
        for j=1:nCha
            for l = 1:nSlices
                x_wavedec_2D{1,j} = wavedec2(real(x{1,j}),waveletStages,waveletFilterName_l12); % atm only based on l1-daubechie
                x_wavedec_2D{1,j+nCha} = wavedec2(imag(x{1,j}),waveletStages,waveletFilterName_l12); % atm only based on l1-daubechie
                x_wave_fine_scale_2D_real = size(x_wavedec_2D{1,j},2) - (3*waveS_l12(waveletStages+1,1)*waveS_l12(waveletStages+1,2));
                threshold_2D(j) = mad(x_wavedec_2D{1,j}(x_wave_fine_scale_2D_real:end),1);
                x_wave_fine_scale_2D_imag = size(x_wavedec_2D{1,j},2) - (3*waveS_l12(waveletStages+1,1)*waveS_l12(waveletStages+1,2));
                threshold_2D(j+nCha) = mad(x_wavedec_2D{1,j}(x_wave_fine_scale_2D_imag:end),1);
            end;
%             x_wavedec_3D{1,j} = wavedec3(real(z{1,j}),waveletStages,waveletFilterName_l1); % atm only based on l1-daubechie
%             x_wavedec_3D{1,j} = wavedec3(imag(z{1,j}),waveletStages,waveletFilterName_l1); % atm only based on l1-daubechie
%             x_wave_fine_scale_3D_real = zeros(0,0);
%             x_wave_fine_scale_3D_imag = zeros(0,0);
%             for k = (7*(x_wavedec_3D{1,j}.level-1)+2):maxWavecells
%                 x_wave_fine_scale_3D_real = [x_wave_fine_scale_3D_real  x_wavedec_3D{1,j}.dec{k+2}];
%                 x_wave_fine_scale_3D_imag = [x_wave_fine_scale_3D_imag  x_wavedec_3D{1,j}.dec{k+2}];
%             end;
%             [size1, size2, size3] = size(x_wave_fine_scale_3D_real);
%             threshold_3D(j) = mad(reshape(x_wave_fine_scale_3D_real,size1*size2*size3,1));
%             threshold_3D(j+nCha) = mad(reshape(x_wave_fine_scale_3D_imag,size1*size2*size3,1));
        end;
    end;
    clear x_wavedec
else
    threshold_2D(1:nCha) = 1;
    threshold_3D(1:nCha) = 1;
end;

thresholdgamma = 0;
for j = 1:nCha
    thresholdgamma = thresholdgamma + (lambdaWave*threshold_2D(j)).^2; % rss of MAD.
end;
thresholdgamma = sqrt(thresholdgamma);

%% initialize metrics:
itr = 0;
metrics.xtime(itr+1)= 0;
x_test = x;
[metrics, ssim_map{1,1}] = get_metrics_itr( im_ref, im_ref_full, At_v, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma, nSlices );

%% recon
for itr = 1:maxitr  % total iter counter

    % landweber step
    for j = 1:nCha
        if strcmp(obj.reconDIM,'3D')
            AtA_x{1,j} = iFFT3D(FFT3D_mask(x_accel{1,j},mask));
        elseif strcmp(obj.reconDIM,'2D')
            AtA_x{1,j} = iFFT2D(FFT2D_mask(x_accel{1,j},mask));
        else
            error('ReconDim is not supported');
        end;
    end;
    
    for j = 1:nCha
        z_coeff{1,j} = ST(x{1,j} - 1/L * (AtA_x{1,j} - At_v{1,j}));  
    end;
    z_coeff_thresh = groupthresh(z_coeff,thresholdgamma,nCha,regularizer);
    for j = 1:nCha
        x{1,j} = BT(z_coeff_thresh{1,j});
    end;
    
    % acceleration
    if flag_fast
        t_old = t_new;
        t_new = (1+sqrt(1+4*t_old^2))/2;   
        for j = 1:nCha        
            x_accel{1,j}=x{1,j}+ (t_old-1)/t_new.*(x{1,j}-x_old{1,j});
        end;
        x_old = x;
    else
        x_accel = x;
    end; 
    
    %% metrics of current itr:
    disp(itr);   
    metrics.xtime(itr+1)= toc(timer_proxA);
    x_test = x_accel;
    [metrics, ssim_map{1,2}] = get_metrics_itr( im_ref, im_ref_full, x_test, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma, nSlices );
    
end;

imageCha = x_test;
for j = 1:nCha
    imageCha{1,j} = turn_image( imageCha{1,j} );
end;
for j = 1:nCha+1
    ssim_map{1,1}{1,j} = turn_image( ssim_map{1,1}{1,j} );
    ssim_map{1,2}{1,j} = turn_image( ssim_map{1,2}{1,j} );
    figure
    imshow(1- max(0, ssim_map{1,1}{1,j}).^4);
    figure
    imshow(1- max(0, ssim_map{1,2}{1,j}).^4);
end;

end