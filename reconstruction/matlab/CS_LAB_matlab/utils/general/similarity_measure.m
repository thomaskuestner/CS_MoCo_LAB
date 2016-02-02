function [ varargout ] = similarity_measure( pic_a, pic_b, varargin )
% measure similarity (intensity-based) between images
% input arguments: 
% pic_a, pic_b      picture a and picture b (2D and 3D); must be the same size!
% grayscale         dynamic range of pixel intensities
% varargin          input similarity measures which should be calculated
%                   output order is the same as input order
%                   'grayscale': grayscale value
%                   'scaleTo': 'own' | 'pic_a' | 'pic_b' scale to max of pic_X
%                   'all': all similarity measures
%                   'basic': smi, snmi, cc, spr, rmse, psnr, ssd, msd, nssd
%                   or: see output arguments
%           
% 
% output arguments: 
% one cell containing all requested similarity measures (in input order)
% or single values:
% smi               Shannon mutual information
% snmi              Shannon normalized mutual information
% efinfo            exclusive F-information
% rmi               Renyi mutual information (additional input ..., 'r', rvalue)
% rnmi              Renyi normalized mutual information (additional input ..., 'r', rvalue)
% tmi               Tsallis mutual information (additional input ..., 'q', qvalue)
% tnmi              Tsallis normalized mutual information (additional input ..., 'q', qvalue)
% ejp               energy of joint probability
% gre               gradient entropy
% finfo             F-information measures (additional input ..., 'alpha', alphavalue)
% ssim/mssim        Structural similarity/Mean SSIM (additional input ..., 'K', [luminance correction, contrast correction], 'window', fspecial(windowType, size, std))
% cc                (Pearson) cross correlation 
% zcc               zero-mean normalized (Pearson) cross correlation coefficient
% spr               Spearman rank correlation
% ket               Kendall's tau
% mse               mean squared error
% nmse              normalized mean squared error
% rmse              root mean squared error
% nrmse             normalized root mean squared error
% psnr              peak signal to noise ratio
% tam               Tanimoto measure
% zca               zero crossings (absolute)
% zcr               zero crossings (relative)
% mr                minimum ratio
% ssd               sum of squared differences
% msd               median of squared differences
% nssd              normalized sum of squared differences
% sad               sum of absolute differences
% zsad              zero-mean sum of absolute differences
% lsad              locally scaled sum of absolute differences
% mad               median of absolute differences
% shd               sum of hamming distance
% besov             besov norm

% all outputs:
% [smi, snmi, efinfo, rmi, rnmi, tmi, tnmi, ejp, gre, finfo, ssim, mssim, cc, zcc, spr, ket, mse, nmse, rmse, nrmse, psnr, tam, zca, zcr, mr, ssd, msd, nssd, sad, zsad, lsad, mad, shd, besov] 

% (c) 2013-2015, Thomas Kuestner, University of Tuebingen,
%     thomas.kuestner@med.uni-tubingen.de


%% set similarity measures
sim_measure_names = {'smi', 'snmi', 'efinfo', 'rmi', 'rnmi', 'tmi', 'tnmi', 'ejp', 'gre', 'finfo', 'ssim', 'mssim', 'cc', 'zcc', 'spr', 'ket', 'mse', 'nmse', 'rmse', 'nrmse', 'psnr', 'tam', 'zca', 'zcr', 'mr', 'ssd', 'msd', 'nssd', 'sad', 'zsad', 'lsad', 'mad', 'shd', 'besov'};
sim_measures = zeros(1,length(sim_measure_names));
if(nargin < 3)
    % basic: smi, snmi, mssim, cc, spr, rmse, psnr, ssd, msd, nssd
    sim_measures([1 2 11 12 14 17 18 23 24 25]) = 1:10;
else
    i = 1;
    input_idx = 1;
    while( i <= nargin-2)
        if(strcmp(varargin{i},'grayscale'))
            grayscale = varargin{i+1};
            i = i+2;
        elseif(strcmp(varargin{i}, 'r'))
            r = varargin{i+1};
            i = i+2;
        elseif(strcmp(varargin{i}, 'q'))
            q = varargin{i+1};
            i = i+2;
        elseif(strcmp(varargin{i}, 'alpha'))
            alpha = varargin{i+1};
            i = i+2;
        elseif(strcmp(varargin{i}, 'K'))
            K = varargin{i+1};
            i = i+2;
        elseif(strcmp(varargin{i}, 'window'))
            window = varargin{i+1};
            i = i+2;
        elseif(strcmp(varargin{i}, 'scaleTo'))
            scaleTo = varargin{i+1};
            i = i+2;
        elseif(strcmp(varargin{i}, 'all'))
            sim_measures = 1:length(sim_measure_names);
            i = i+1;
        elseif(strcmp(varargin{i}, 'basic'))
            sim_measures([1 2 11 12 14 17 18 23 24 25]) = input_idx:input_idx+9;
            i = i+1;
        else
            sim_measures(strcmpi(varargin{i},sim_measure_names)) = input_idx;
            input_idx = input_idx + 1;
            i = i+1;
        end            
    end
end

    
%% set #bins
if(~exist('grayscale','var') || isempty(grayscale))
    grayscale = 256;
end
x=0:grayscale-1; 
    

%% images
if(any(size(pic_a) ~= size(pic_b)))
    error('Similarity_measure: Input images have different sizes');
end

if(~exist('scaleTo','var') || isempty(scaleTo))
    scaleTo = 'own';
end
% first scale dynamic ranges of images to fit in {0,...,grayscale-1}
% otherwise there will be an inconsistency of the results
if(strcmp(scaleTo,'own'))
    pic_a=(pic_a-min(pic_a(:)))*(x(end)-x(1))./(max(pic_a(:))-min(pic_a(:))); 
    pic_b=(pic_b-min(pic_b(:)))*(x(end)-x(1))./(max(pic_b(:))-min(pic_b(:))); 
elseif(strcmp(scaleTo,'pic_a'))
    maxScale = max(pic_a(:));
    pic_a=(pic_a-min(pic_a(:)))*(x(end)-x(1))./(maxScale-min(pic_a(:))); 
    pic_b=(pic_b-min(pic_b(:)))*(x(end)-x(1))./(maxScale-min(pic_b(:)));
elseif(strcmp(scaleTo,'pic_b'))
    maxScale = max(pic_b(:));
    pic_a=(pic_a-min(pic_a(:)))*(x(end)-x(1))./(maxScale-min(pic_a(:))); 
    pic_b=(pic_b-min(pic_b(:)))*(x(end)-x(1))./(maxScale-min(pic_b(:)));
end


%% marginal histograms
n_a = histc(pic_a(:),x);
n_b = histc(pic_b(:),x);

% probabilities
p_a = n_a ./ sum(n_a);
p_b = n_b ./ sum(n_b);

%% joint histogram
n_ab = jointHistogram(pic_a, pic_b, grayscale, 1);

%% joint probability
p_ab = n_ab ./ sum(sum(n_ab));

%% Shannon mutual informations (SMI/SNMI) 
if(sim_measures(strcmp('smi',sim_measure_names)) || sim_measures(strcmp('snmi',sim_measure_names)))
    % marginal entropies
    H_a = -sum(p_a(p_a ~= 0) .* log2(p_a(p_a ~= 0)));
    H_b = -sum(p_b(p_b ~= 0) .* log2(p_b(p_b ~= 0)));
    % joint entropy
    H_ab = -sum(p_ab(p_ab ~= 0) .* log2(p_ab(p_ab ~= 0)));
end

%% Shannon mutual information (SMI)
if(sim_measures(strcmp('smi',sim_measure_names)))
    smi = H_a + H_b - H_ab;
end

%% Shannon normalized mutual information (SNMI)
if(sim_measures(strcmp('snmi',sim_measure_names)))
    snmi = (H_a + H_b)/(H_ab);
end


%% exclusive F-information (EFINFO)
if(sim_measures(strcmp('efinfo',sim_measure_names)))
    efinfo = 2 * H_ab - H_a - H_b;
end


%% Renyi mutual informations (RMI/RNMI)
if(sim_measures(strcmp('rmi',sim_measure_names)) || sim_measures(strcmp('rnmi',sim_measure_names)))
    % needs variable r
    if(~exist('r','var') || isempty(r))
        r = 0.5;
    end
    formula = '1/(1-r) .* log2(sum((p_X(p_X ~= 0)).^r)/sum(p_X(p_X ~= 0)));';
    for idx = {'a','b','ab'}
        eval(sprintf('H_%s = %s',idx{1},strrep(formula,'X',idx{1})));
    end
    if(sim_measures(strcmp('rmi',sim_measure_names)))
        rmi = H_a + H_b - H_ab;
    end
    if(sim_measures(strcmp('rnmi',sim_measure_names)))
        rnmi = (H_a + H_b)/(H_ab);
    end
end


%% Tsallis mutual informations (TMI/TNMI)
if(sim_measures(strcmp('tmi',sim_measure_names)) || sim_measures(strcmp('tnmi',sim_measure_names)))
    % needs variable q
    if(~exist('q','var') || isempty(q))
        q = 0.5;
    end
    formula = '1/(q-1) .* (1 - sum((p_X(p_X ~= 0)).^r));';
    for idx = {'a','b','ab'}
        eval(sprintf('H_%s = %s',idx{1},strrep(formula,'X',idx{1})));
    end
    if(sim_measures(strcmp('tmi',sim_measure_names)))
        tmi = H_a + H_b + (1-q)*H_a*H_b - H_ab;
    end
    if(sim_measures(strcmp('tnmi',sim_measure_names)))
        tnmi = (H_a + H_b + (1-q)*H_a*H_b)/(H_ab);
    end
end

%% energy of joint probability (EJP)
if(sim_measures(strcmp('ejp',sim_measure_names)))
    %jpd = n_ab./(numel(pic_a(:))); % joint probability distribution, same as p_ab
    jpd = p_ab;
    ejp = sum(sum(sum((jpd).^2)));
end


%% Gradient entropy (GRE)
if(sim_measures(strcmp('gre',sim_measure_names)))
%     sobel(:,:,1) = [1 2 1; 2 4 2; 1 2 1];
%     sobel(:,:,2) = zeros(3,3);
%     sobel(:,:,3) = -sobel(:,:,1);
%     G_a = abs(convn(pic_a,sobel),'same');
%     G_b = abs(convn(pic_b,sobel),'same');
    
    G_a = zeros(size(pic_a));
    G_b = G_a;
    G_ab = G_a;
    for i=1:size(pic_a,3)
        G_a(:,:,i) = abs(edge(pic_a(:,:,i),'sobel'));
        G_b(:,:,i) = abs(edge(pic_b(:,:,i),'sobel'));
        G_ab(:,:,i) = abs(edge(abs(pic_a(:,:,i) - pic_b(:,:,i)),'sobel'));
    end

%     G_ab = abs(G_a - G_b);
    G_a = G_a./sum(G_a(:));
    G_b = G_b./sum(G_b(:));
    G_ab = G_ab./sum(G_ab(:));
    
    H_a = -sum(G_a(G_a ~= 0) .* log2(G_a(G_a ~= 0)));
    H_b = -sum(G_b(G_b ~= 0) .* log2(G_b(G_b ~= 0)));
    H_ab = -sum(G_ab(G_ab ~= 0) .* log2(G_ab(G_ab ~= 0)));
%     gre = H_ab;
    
    gre = 0.5*(H_a + H_b)/H_ab;
%     gre = abs(H_a - H_b);
%     gre = H_a + H_b - H_ab;

end


%% F-information measures
if(sim_measures(strcmp('finfo',sim_measure_names)))
    % needs variable alpha
    if(~exist('alpha','var') || isempty(alpha))
        alpha = 0.5;
    end
    if(alpha ~= 0 && alpha ~= 1)
        finfo.Ia = 1/(alpha * (alpha -1)) .* sum(sum((p_ab).^alpha ./(p_a * p_b').^(alpha-1) - 1));
    else
        finfo.Ia = NaN;
    end
    if((0 < alpha) && (alpha <= 1))
        finfo.Ma = sum(sum((abs((p_ab).^alpha - (p_a * p_b').^alpha)).^(1/alpha)));
    else
        finfo.Ma = NaN;
    end
    if(alpha > 1)
        finfo.Xa = sum(sum(((abs(p_ab - p_a * p_b')).^alpha)./((p_a * p_b').^(alpha-1))));
    else
        finfo.Xa = NaN;
    end 
end


%% Structural similarity (SSIM)
if(sim_measures(strcmp('ssim',sim_measure_names)) || sim_measures(strcmp('mssim',sim_measure_names)))
    if(~exist('K','var') || isempty(K))
        K = [0.01 0.03];
    end
    if(~exist('window','var') || isempty(window))
        window = fspecial('gaussian', 11, 1.5);
    end
    if(ndims(pic_a) <= 3)
        mssim_tmp = zeros(1,size(pic_a,3));
        ssim = zeros(size(pic_a,1) - size(window,1) + 1, size(pic_a,2) - size(window,2) + 1, size(pic_a,3));
        for i=1:size(pic_a,3)
            [mssim_tmp(i), ssim(:,:,i)] = ssim_index(pic_a(:,:,i), pic_b(:,:,i), K, window, max(x));
        end
        mssim = mean(mssim_tmp);
    end
end


%% (Pearson) normalized cross correlation coefficient (NCC/CC)
% Hint: Matlab uses default 1/(n-1) for expectation value instead of 1/n,
% but in correlation coefficient this cancels out

% for 1/(n-1) and 1/n
if(sim_measures(strcmp('cc',sim_measure_names)))
    cc = corrcoef(pic_a,pic_b);
    cc = cc(1,2);
end

%% zero-mean normalized cross correlation coefficient (ZCC)
if(sim_measures(strcmp('zcc',sim_measure_names)))
    zcc = sum(sum(sum(pic_a.*pic_b)))/sqrt(sum(sum(sum((pic_a).^2))) .* sum(sum(sum((pic_b).^2))));
end


%% Spearman rank correlation (SPR)
if(sim_measures(strcmp('spr',sim_measure_names)))
    spr = corr(pic_a(:), pic_b(:), 'type', 'Spearman');
%     spr = 1 - ((6*sum((pic_a(:) - pic_b(:)).^2))/(numel(pic_a(:)) * (numel(pic_a(:)).^2 - 1)));
end


%% Kendall's tau (KET)
if(sim_measures(strcmp('ket',sim_measure_names)))
    ket = corr(pic_a(:), pic_b(:), 'type', 'Kendall');
end


%% mean squared error (MSE)
if(sim_measures(strcmp('mse',sim_measure_names)))
    mse = sum((pic_a(:) - pic_b(:)).^2)/numel(pic_a(:));
end


%% normalized mean squared error (NMSE)
if(sim_measures(strcmp('nmse',sim_measure_names)))
    nmse = sum((pic_a(:) - pic_b(:)).^2)/(sum(pic_b(:).^2));
end


%% root mean squared error (RMSE)
if(sim_measures(strcmp('rmse',sim_measure_names)))
    rmse = sqrt(sum(abs(pic_a(:) - pic_b(:)).^2)/numel(pic_a(:)));
end


%% normalized root mean squared error (NRMSE)
if(sim_measures(strcmp('nrmse',sim_measure_names)))
    nrmse = 1/(max(abs(pic_a(:)))-min(abs(pic_a(:)))) .* sqrt(1/numel(pic_a(:)) .* sum(abs(pic_a(:) - pic_b(:)).^2));
end


%% peak to signal noise ratio (PSNR)
if(sim_measures(strcmp('psnr',sim_measure_names)))
    if(exist('mse','var'))
        psnr = 10*log10(grayscale^2 / abs(mse));
    else
        psnr = 10*log10(grayscale^2 / abs(sum(pic_a(:) - pic_b(:))/numel(pic_a(:))));
    end
end


%% Tanimoto measure (TAM)
if(sim_measures(strcmp('tam',sim_measure_names)))
    tam = (pic_a(:)' * pic_b(:))/(sum(abs(pic_a(:) - pic_b(:)).^2) + pic_a(:)' * pic_b(:));
end


%% zero crossings (absolute and relative) (ZCA/ZCR)
if(sim_measures(strcmp('zca',sim_measure_names)) || sim_measures(strcmp('zcr',sim_measure_names)))
    pic_diff = pic_a - pic_b;
    pic_diff = pic_diff(:);
    if(sim_measures(strcmp('zca',sim_measure_names)))
        zca = length(find(pic_diff(1:end-1).*pic_diff(2:end) < 0));
    end
    if(sim_measures(strcmp('zcr',sim_measure_names)))
        zcr = zca/(numel(pic_diff(:)));
    end
end


%% minimum ratio (MR)
if(sim_measures(strcmp('mr',sim_measure_names)))
    mr = sum(sum(sum(min(pic_a,pic_b))))/numel(pic_a(:));
end
    

%% sum of squared differences (SSD)
if(sim_measures(strcmp('ssd',sim_measure_names)))
    ssd = sqrt(sum(sum(sum((pic_a-pic_b).^2))))/(numel(pic_a));
end


%% median of squared differences (MSD)
if(sim_measures(strcmp('msd',sim_measure_names)))
    msd = median(median(median((pic_a - pic_b).^2)))/(numel(pic_a));
end


%% normalized sum of squared differences (nssd)
if(sim_measures(strcmp('nssd',sim_measure_names)))
    nssd = sqrt(sum(sum(sum(((pic_a - mean(mean(mean(pic_a))))/std(std(std(pic_a,1),1),1) - (pic_b - mean(mean(mean(pic_b))))/std(std(std(pic_b,1),1),1)).^2))));
end


%% sum of absolute differences (SAD)
if(sim_measures(strcmp('sad',sim_measure_names)))
    sad = sum(sum(sum(abs(pic_a - mean(pic_a(:)) - pic_b + mean(pic_b(:))))))/(numel(pic_a));
end
    

%% zero-mean sum of absolute differences (ZSAD)
if(sim_measures(strcmp('zsad',sim_measure_names)))
    zsad = sum(sum(sum(abs(pic_a - pic_b))))/(numel(pic_a));
end
    

%% locally scaled sum of absolute differences (LSAD)
if(sim_measures(strcmp('lsad',sim_measure_names)))
    lsad = sum(sum(sum(abs(pic_a - (mean(pic_a(:))/mean(pic_b(:)) * pic_b)))))/(numel(pic_a));
end

    
%% median of absolute differences (MAD)
if(sim_measures(strcmp('mad',sim_measure_names)))
    mad = median(median(median(pic_a - pic_b)))/(numel(pic_a));
end


%% sum of hamming distance (SHD)
if(sim_measures(strcmp('shd',sim_measure_names)))
    shd = sum(pic_a(:) == pic_b(:))/(numel(pic_a));
end


%% besov norm (BESOV)
if(sim_measures(strcmp('besov',sim_measure_names)))
    besov = besovNorm( pic_a - pic_b );
end

%% outputs
if(nargout == 0)
    % just print results
    for i=1:nnz(sim_measures)
        idx = sim_measures == min(sim_measures(sim_measures ~= 0)); % for given input order
        if(exist(sim_measure_names{idx},'var'))
            if(strcmp(sim_measure_names{idx},'finfo'))
                fprintf('F.Ia: \t%.4f\nF.Ma: \t%.4f\nF.Xa: \t%.4f\n', finfo.Ia, finfo.Ma, finfo.Xa);
            elseif(strcmp(sim_measure_names{idx},'ssim'))
                fprintf('SSIM: Check Map\n');
            else
                fprintf('%s: \t%.4f\n', upper(sim_measure_names{idx}), eval(sim_measure_names{idx}));
            end
        end
        sim_measures(idx) = 0;
    end
elseif(nargout == nnz(sim_measures))
    for i=1:nargout
        idx = sim_measures == min(sim_measures(sim_measures ~= 0)); % for given input order
        if(exist(sim_measure_names{idx},'var'))
            eval(sprintf('varargout{i} = %s;', sim_measure_names{idx}));
        end
        sim_measures(idx) = 0;
    end    
else
    outMeasures = cell(1,nnz(sim_measures));
    for i=1:nnz(sim_measures)
        idx = sim_measures == min(sim_measures(sim_measures ~= 0)); % for given input order
        if(exist(sim_measure_names{idx},'var'))
            eval(sprintf('outMeasures{i} = %s;', sim_measure_names{idx}));
        end
        sim_measures(idx) = 0;
    end
    varargout{1} = outMeasures;
end
end


function n_ab = jointHistogram(pic_a, pic_b,grayscale,ver)
% calculate a joint histogram

if(ver == 1)
    xi = linspace(min(pic_a(:)),max(pic_a(:)),grayscale);
    yi = linspace(min(pic_b(:)),max(pic_b(:)),grayscale);

    xr = interp1(xi,1:numel(xi),pic_a(:),'nearest');
    yr = interp1(yi,1:numel(yi),pic_b(:),'nearest');

    n_ab = accumarray([xr yr], 1, [grayscale grayscale]);

else
    x=0:grayscale-1; 
    n_ab=zeros(grayscale); 
    for id=0:grayscale-1 
        n_ab(id+1,:) = histc(pic_b(pic_a==id),x,1); 
    end
end
end

function [mssim, ssim_map] = ssim_index(img1, img2, K, window, L)

%========================================================================
%SSIM Index, Version 1.0
%Copyright(c) 2003 Zhou Wang
%All Rights Reserved.
%
%The author is with Howard Hughes Medical Institute, and Laboratory
%for Computational Vision at Center for Neural Science and Courant
%Institute of Mathematical Sciences, New York University.
%
%----------------------------------------------------------------------
%Permission to use, copy, or modify this software and its documentation
%for educational and research purposes only and without fee is hereby
%granted, provided that this copyright notice and the original authors'
%names appear on all copies and supporting documentation. This program
%shall not be used, rewritten, or adapted as the basis of a commercial
%software or hardware product without first obtaining permission of the
%authors. The authors make no representations about the suitability of
%this software for any purpose. It is provided "as is" without express
%or implied warranty.
%----------------------------------------------------------------------
%
%This is an implementation of the algorithm for calculating the
%Structural SIMilarity (SSIM) index between two images. Please refer
%to the following paper:
%
%Z. Wang, A. C. Bovik, H. R. Sheikh, and E. P. Simoncelli, "Image
%quality assessment: From error measurement to structural similarity"
%IEEE Transactios on Image Processing, vol. 13, no. 1, Jan. 2004.
%
%Kindly report any suggestions or corrections to zhouwang@ieee.org
%
%----------------------------------------------------------------------
%
%Input : (1) img1: the first image being compared
%        (2) img2: the second image being compared
%        (3) K: constants in the SSIM index formula (see the above
%            reference). defualt value: K = [0.01 0.03]
%        (4) window: local window for statistics (see the above
%            reference). default widnow is Gaussian given by
%            window = fspecial('gaussian', 11, 1.5);
%        (5) L: dynamic range of the images. default: L = 255
%
%Output: (1) mssim: the mean SSIM index value between 2 images.
%            If one of the images being compared is regarded as 
%            perfect quality, then mssim can be considered as the
%            quality measure of the other image.
%            If img1 = img2, then mssim = 1.
%        (2) ssim_map: the SSIM index map of the test image. The map
%            has a smaller size than the input images. The actual size:
%            size(img1) - size(window) + 1.
%
%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   [mssim ssim_map] = ssim_index(img1, img2);
%
%Advanced Usage:
%   User defined parameters. For example
%
%   K = [0.05 0.05];
%   window = ones(8);
%   L = 100;
%   [mssim ssim_map] = ssim_index(img1, img2, K, window, L);
%
%See the results:
%
%   mssim                        %Gives the mssim value
%   imshow(max(0, ssim_map).^4)  %Shows the SSIM index map
%
%========================================================================


if (nargin < 2 || nargin > 5)
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(img1) ~= size(img2))
   ssim_index = -Inf;
   ssim_map = -Inf;
   return;
end

[M N] = size(img1);

if (nargin == 2)
   if ((M < 11) || (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);	%
   K(1) = 0.01;								      % default settings
   K(2) = 0.03;								      %
   L = 255;                                  %
end

if (nargin == 3)
   if ((M < 11) || (N < 11))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 4)
   [H W] = size(window);
   if ((H*W) < 4 || (H > M) || (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   L = 255;
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 5)
   [H W] = size(window);
   if ((H*W) < 4 || (H > M) || (W > N))
	   ssim_index = -Inf;
	   ssim_map = -Inf;
      return
   end
   if (length(K) == 2)
      if (K(1) < 0 || K(2) < 0)
		   ssim_index = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   ssim_index = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));
img1 = double(img1);
img2 = double(img2);

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
	denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

end


function out = besovNorm( img, p, beta, gamma, filter, level )
%BESOVNORM calculate besov norm

if(nargin < 2)
    p = 1;
    beta = 4/9;
    gamma = 4;
    filter = 'db1';
    level = 3;
elseif(nargin < 5)
    filter = 'db1';
    level = 3;
end
    

if(ndims(img) == 3)
    imgWave = wavedec3(img,level,filter,'mode','zpd');
    imgWave = imgWave.dec;
    step = 7;
elseif(ismatrix(img))
    imgWave = wavedec2(img,level,filter,'mode','zpd');
    step = 4;
else
    error('besovNorm(): Unknown matrix dimensions');
end
    

position = length(imgWave);
out = 0;
for j=1:level
    if(ndims(img) == 3)
        tmp = abs(cell2mat(shiftdim(imgWave(position:-1:position-step+1,:),-2))).^beta;
    else
        tmp = abs(imgWave(position:-1:position-step+1)).^beta;
    end
    tmp = (sum(tmp(:)))^(1/beta);
    out = out + (2^(-j * (p + 1/2 - 1/beta)) * tmp)^gamma;
    position = position - step;
end
% approximation of level N
if(ndims(img) == 3)
    tmp = abs(cell2mat(shiftdim(imgWave(1,1),-2))).^beta;
else
    tmp = abs(imgWave(1)).^beta;
end
tmp = (sum(tmp(:)))^(1/beta);
out = out + (2^(-(level+1) * (p + 1/2 - 1/beta)) * tmp)^gamma;
out = out^(1/gamma);

end


