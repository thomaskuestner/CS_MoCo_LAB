function [imageCha, ssim_map, metrics] = algo_A_PG( obj,input,TransformSpecifics )
%%
% Accelerated Proximal Gradient with variable regularizer

% M. Fischer, April 2016

%%
ssim_map = [];
metrics = [];
timer_proxA = tic;

%% variables:
if strcmp(obj.reconDIM,'3D')
    if strcmp(TransformSpecifics.TransformName,'dwt') || strcmp(TransformSpecifics.TransformName,'cplxdt') || strcmp(TransformSpecifics.TransformName,'shearlet')
        flag_cells = true;
    else
        flag_cells = false;
    end;
elseif strcmp(obj.reconDIM,'2D')
    if strcmp(TransformSpecifics.TransformName,'cplxdt') || strcmp(TransformSpecifics.TransformName,'cplxdt_matlab') || strcmp(TransformSpecifics.TransformName,'shearlet')
        flag_cells = true;
    else
        flag_cells = false;
    end;
end;

% from obj:
% maxitr = obj.maxitr;
maxitr = obj.iNINNER;
% n1 = obj.measPara.dim(1);
% n2 = obj.measPara.dim(2);
n1 = input.n1;
n2 = input.n2;
if strcmp(obj.reconDIM,'3D')
    n3 = obj.measPara.dim(3);
elseif strcmp(obj.reconDIM,'2D')
    n3 = 1;
end;
% nSlices = obj.measPara.dim(3);
nCha = obj.measPara.dim(5);
lambdaWave = obj.lambda;

L = 2;
t_new = 1;
if strcmp(input.regularizer,'l1')
    flag_fast = 1;
    L = 1.5; % pFISTA is an approximate scheme. Increase L to ensure empirical convergence.
else
    flag_fast = 0;
    L = 1;
end;

%% operators
ST = @(x) operator_frame(x,n1,n2,n3,1,0,0,...
    @(x) sparse_trans(x,TransformSpecifics,obj.reconDIM));
BT = @(x) operator_frame(x,n1,n2,n3,1,0,1,...
    @(x) back_trans(x,TransformSpecifics,obj.reconDIM,n1,n2,n3));


% from input:
b=input.b;
mask = input.mask;
regularizer = input.regularizer;

% im_ref = input.im_ref;
% if strcmp(obj.reconDIM,'3D')
%     im_ref_full = zeros(n1,n2,n3);
% elseif strcmp(obj.reconDIM,'2D')
%     im_ref_full = zeros(n1,n2);
% end;
% 
% for j = 1:nCha
%     im_ref_full = im_ref_full + abs(im_ref{1,j}).^2;
% end;
% im_ref_full = sqrt(im_ref_full);
clear input

%% initialize cells/vectors
for j=1:nCha
    if strcmp(obj.reconDIM,'3D')
        A_trans_v{1,j} = iFFT3D(b{1,j});
    elseif strcmp(obj.reconDIM,'2D')
        A_trans_v{1,j} = iFFT2D(b{1,j});
    else
        error('ReconDim is not supported');
    end;
end;

% starting point:
x = A_trans_v; % (x0 = A_trans_v)
x_accel = x;
x_old = x;

AtA_x =  cell(1,nCha);
for j = 1:nCha
        AtA_x{1,j} = 0;
end;

if ~flag_cells
    z_coeff = AtA_x;  
    z_coeff_thresh = AtA_x;
else
    z_coeff = cell(TransformSpecifics.CellAmount,1);
    for k = 1:TransformSpecifics.CellAmount
        z_coeff{k} = cell(1,nCha);
        for j = 1:nCha
            z_coeff{k}{1,j} = 0;
        end;
    end;
    z_coeff_thresh = z_coeff;
    
    z_coeff_helper = cell(1,nCha);
    for j = 1:nCha
        z_coeff_helper{1,j} = cell(TransformSpecifics.CellAmount,1);
        for k = 1:TransformSpecifics.CellAmount
            z_coeff_helper{1,j}{k} = 0;
        end;
    end;
    z_coeff_thresh_helper = z_coeff_helper;
end;

%% MAD dependent lambda:
flag_MAD = true;
if flag_MAD
    threshold = get_MAD(x,nCha,TransformSpecifics);
else
    threshold(1:nCha) = 1;
end;

thresholdgamma = 0;
for j = 1:nCha
    thresholdgamma = thresholdgamma + ((lambdaWave*threshold(j))).^2; % rss of MAD.
end;
thresholdgamma = sqrt(thresholdgamma)*L;

%% initialize metrics:
    itr = 0;
%    	metrics.xtime(itr+1)= 0;    
%     [metrics, ssim_map{1,1}] = get_metrics_itr( im_ref, im_ref_full, x, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma, n3 );

%% recon
dispProgress('Proximal Average', 0, maxitr);
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
    
    if ~flag_cells
        for j = 1:nCha
            z_coeff{1,j} = ST(x{1,j} - 1/L * (AtA_x{1,j} - A_trans_v{1,j}));
        end;
        z_coeff_thresh = groupthresh(z_coeff,thresholdgamma,nCha,regularizer);
        for j = 1:nCha
            x{1,j} = BT(z_coeff_thresh{1,j});
        end;
    else
        for j = 1:nCha
            z_coeff_helper{1,j} = ST(x{1,j} - 1/L * (AtA_x{1,j} - A_trans_v{1,j}));
            for k = 1:TransformSpecifics.CellAmount% TransformSpecifics.Stages+1
                z_coeff{k}{1,j} = z_coeff_helper{1,j}{k}; % switch for groupthresh
            end;
        end;
        for k = 1:TransformSpecifics.CellAmount
            z_coeff_thresh{k} = groupthresh(z_coeff{k},thresholdgamma,nCha,regularizer);
        end;
        for j = 1:nCha
            for k = 1:TransformSpecifics.CellAmount
                z_coeff_thresh_helper{1,j}{k} = z_coeff_thresh{k}{1,j}; % switch for backtrafo
            end;
            x{1,j} = BT(z_coeff_thresh_helper{1,j});
        end;
    end;
    
    % acceleration
    if flag_fast
        t_old = t_new;
        t_new = (1+sqrt(1+4*t_old^2))/2;
        for j = 1:nCha
            x_accel{1,j} = x{1,j} + ((t_old-1)/t_new)*(x{1,j} - x_old{1,j});
        end;
        x_old = x;
    else
        x_accel = x;
    end;
    
    dispProgress('Proximal Average', itr/maxitr);
    
    %% metrics of current itr:
%     disp(itr);
    
%     metrics.xtime(itr+1)= toc(timer_proxA);
%     [metrics, ssim_map{1,2}] = get_metrics_itr( im_ref, im_ref_full, x, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma, n3 );
end;
dispProgress('Proximal Average', 'Close');

imageCha = x;
% for j = 1:nCha
%     imageCha{1,j} = turn_image( imageCha{1,j} );
% end;
% for j = 1:nCha+1
%     ssim_map{1,1}{1,j} = turn_image( ssim_map{1,1}{1,j} );
%     ssim_map{1,2}{1,j} = turn_image( ssim_map{1,2}{1,j} );
%     figure
%     imshow(1- max(0, ssim_map{1,1}{1,j}(:,:,ceil(n3/2))).^4);
%     figure
%     imshow(1- max(0, ssim_map{1,2}{1,j}(:,:,ceil(n3/2))).^4);
% end;

end