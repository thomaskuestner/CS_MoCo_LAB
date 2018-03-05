function [imageCha, ssim_map, metrics] = algo_A_ADMM( obj,input,TransformSpecifics )  
%%
% ADMM with variable regularizer

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
    n3 = obj.measPara.dim(2);
elseif strcmp(obj.reconDIM,'2D')
    n3 = 1;
end;
% nSlices = obj.measPara.dim(3);
nCha = obj.measPara.dim(5);
lambdaWave = obj.lambda;
tau = obj.mue;

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
rhs = A_trans_v; % different start (x0 = F^-1 K^-1 F A_trans_v); %used atm
% rhs: right hand side

rho_rhs =  cell(1,nCha);
for j = 1:nCha
        rho_rhs{1,j} = 0;
end;
K_inv = rho_rhs;

if ~flag_cells
    z_coeff = rho_rhs;
    d_aux_helper = rho_rhs;
    d_aux = rho_rhs;
    d_aux_old = rho_rhs;
    gamma = rho_rhs;
    gamma_old = rho_rhs;
else
    z_coeff = cell(TransformSpecifics.CellAmount,1);
    for k = 1:TransformSpecifics.CellAmount
        z_coeff{k} = cell(1,nCha);
        for j = 1:nCha
            z_coeff{k}{1,j} = 0;
        end;
    end;
    d_aux_helper = z_coeff;
    d_aux = z_coeff;
    d_aux_old = z_coeff;
    gamma = z_coeff;
    gamma_old = z_coeff;
    
    z_coeff_helper = cell(1,nCha);
    for j = 1:nCha
        z_coeff_helper{1,j} = cell(TransformSpecifics.CellAmount,1);
        for k = 1:TransformSpecifics.CellAmount
            z_coeff_helper{1,j}{k} = 0;
        end;
    end;
    bt_helper = z_coeff_helper;
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
    thresholdgamma = thresholdgamma + ((lambdaWave*threshold(j))/tau).^2; % rss of MAD.
end;
thresholdgamma = sqrt(thresholdgamma);

%% initialize metrics:
    itr = 0;
%    	metrics.xtime(itr+1)= 0;    
%     [metrics, ssim_map{1,1}] = get_metrics_itr( im_ref, im_ref_full, x, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma, n3 );
    
%% recon
% inverse for rho - step:
for j = 1:nCha
    K_inv{1,j} = 1./(mask + tau.*ones(n1,n2,n3));
end;

dispProgress('ADMM', 0, maxitr);
for itr = 1:maxitr  % total iter counter
    %% rho - step
    if itr > 1 % for i = 1 already set (x{1,j} = A_trans_v{1,j})
        % rhs part:
        if ~flag_cells
            for j = 1:nCha
                rho_rhs{1,j} = BT(d_aux_old{1,j} - gamma_old{1,j}./tau);
            end
        else
            for j = 1:nCha
                for k = 1:TransformSpecifics.CellAmount
                    bt_helper{1,j}{k} = d_aux_old{k}{1,j} - gamma_old{k}{1,j}./tau; % switch for backtrafo
                end;
                rho_rhs{1,j} = BT(bt_helper{1,j});
            end;
        end;
        
        for j = 1:nCha
            rhs{1,j} = A_trans_v{1,j} + tau .* rho_rhs{1,j};
            
            if strcmp(obj.reconDIM,'3D')
                x{1,j} = iFFT3D(K_inv{1,j}.*FFT3D(rhs{1,j}));
            elseif strcmp(obj.reconDIM,'2D')
                x{1,j} = iFFT2D(K_inv{1,j}.*FFT2D(rhs{1,j}));
            else
                error('ReconDim is not supported');
            end;
        end;
    end;
    
    %% auxiliary variable
    if ~flag_cells
        for j = 1:nCha    
            z_coeff{1,j} = ST(x{1,j});
            d_aux_helper{1,j} = z_coeff{1,j} + gamma_old{1,j}./tau;
        end;
        
        d_aux = groupthresh(d_aux_helper,thresholdgamma,nCha,regularizer);
        
        % Lagrange multiplier        
        for j = 1:nCha
           gamma{1,j} = gamma_old{1,j} + tau .* (z_coeff{1,j} - d_aux{1,j});
        end;
    else
        for j = 1:nCha      
            z_coeff_helper{1,j} = ST(x{1,j});
            for k = 1:TransformSpecifics.CellAmount% TransformSpecifics.Stages+1
                z_coeff{k}{1,j} = z_coeff_helper{1,j}{k}; % switch for groupthresh
                d_aux_helper{k}{1,j} = z_coeff{k}{1,j} + gamma_old{k}{1,j}./tau;
            end;
        end;
        
        for k = 1:TransformSpecifics.CellAmount
            d_aux{k} = groupthresh(d_aux_helper{k},thresholdgamma,nCha,regularizer);
        end;
        
        % Lagrange multiplier            
        for k = 1:TransformSpecifics.CellAmount
            for j = 1:nCha
                gamma{k}{1,j} = gamma_old{k}{1,j} + tau .* (z_coeff{k}{1,j} - d_aux{k}{1,j});
            end;
        end;
    end;
    
    d_aux_old = d_aux;
    gamma_old = gamma;
    
    dispProgress('ADMM', itr/maxitr);

    %% metrics of current itr:
%     disp(itr);
%     
%     metrics.xtime(itr+1)= toc(timer_proxA);
%     [metrics, ssim_map{1,2}] = get_metrics_itr( im_ref, im_ref_full, x, itr, maxitr, nCha, n1, n2, metrics, obj.K_1, obj.K_2, obj.W_size, obj.W_sigma, n3 );
end;
dispProgress('ADMM', 'Close');

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

