function [ metrics, ssim_map ] = get_metrics_itr( im_ref, im_ref_full, im_complex, itr, maxitr, nCha, n1, n2, metrics, K_1, K_2, W_size, W_sigma, nSlices )
% gathers metrics
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

if nargin < 14
    nSlices = 1;
end;

flag_plotMap = false;

    for j = 1:nCha+1
        if j <= nCha
            im_y = abs(im_complex{1,j});
            metrics.snr{1,j}(itr+1) = snr(im_y, im_ref{1,j}); % from WaTMRI
            metrics.psnr{1,j}(itr+1) = similarity_measure( im_ref{1,j}, im_y, 'psnr' );
            [metrics.ssim{1,j}(itr+1), ~] = similarity_measure_ssim( im_ref{1,j}, im_y, K_1, K_2, W_size, W_sigma );
            if itr < 1
                [~, ssim_map{1,j}] = similarity_measure_ssim( im_ref{1,j}, im_y, K_1, K_2, W_size, W_sigma );
                ssim_map{1,j} = padarray(ssim_map{1,j},[floor(W_size/2) floor(W_size/2)]);
            elseif itr >= maxitr
                [~, ssim_map{1,j}] = similarity_measure_ssim( im_ref{1,j}, im_y, K_1, K_2, W_size, W_sigma );
                ssim_map{1,j} = padarray(ssim_map{1,j},[floor(W_size/2) floor(W_size/2)]);
            else
                ssim_map = 0;
            end;
            if flag_plotMap
                figure
                imshow(1-max(0, metrics.ssim_map{1,j}{1,1}).^4);
            end;
            metrics.snmi{1,j}(itr+1) = similarity_measure( im_ref{1,j}, im_y, 'snmi' );
            metrics.cnmi{1,j}(itr+1) = similarity_measure( im_ref{1,j}, im_y, 'cnmi' );
            metrics.vi{1,j}(itr+1) = similarity_measure( im_ref{1,j}, im_y, 'vi' );
            metrics.cc{1,j}(itr+1) = similarity_measure( im_ref{1,j}, im_y, 'cc' );
            metrics.nrmse{1,j}(itr+1) = similarity_measure( im_ref{1,j}, im_y, 'nrmse' );
        else
            im_y = zeros(n1,n2,nSlices);
            for h = 1:nCha
                im_y = im_y + abs(im_complex{1,h}).^2;
            end;
            im_y = sqrt(im_y);
            metrics.snr{1,j}(itr+1) = snr(im_y, im_ref_full); % from WaTMRI
            metrics.psnr{1,j}(itr+1) = similarity_measure( im_ref_full, im_y, 'psnr' );
            [metrics.ssim{1,j}(itr+1), ~] = similarity_measure_ssim( im_ref_full, im_y, K_1, K_2, W_size, W_sigma );
            if itr < 1
                [~, ssim_map{1,j}] = similarity_measure_ssim( im_ref_full, im_y, K_1, K_2, W_size, W_sigma );
                ssim_map{1,j} = padarray(ssim_map{1,j},[floor(W_size/2) floor(W_size/2)]);
            elseif itr >= maxitr
                [~, ssim_map{1,j}] = similarity_measure_ssim( im_ref_full, im_y, K_1, K_2, W_size, W_sigma );
                ssim_map{1,j} = padarray(ssim_map{1,j},[floor(W_size/2) floor(W_size/2)]);
            end;
            if flag_plotMap
                figure
                imshow(1-max(0, metrics.ssim_map{1,j}{1,1}).^4);
            end;
            metrics.snmi{1,j}(itr+1) = similarity_measure( im_ref_full, im_y, 'snmi' );
            metrics.cnmi{1,j}(itr+1) = similarity_measure( im_ref_full, im_y, 'cnmi' );
            metrics.vi{1,j}(itr+1) = similarity_measure( im_ref_full, im_y, 'vi' );
            metrics.cc{1,j}(itr+1) = similarity_measure( im_ref_full, im_y, 'cc' );
            metrics.nrmse{1,j}(itr+1) = similarity_measure( im_ref_full, im_y, 'nrmse' );
        end;
    end;
end

