function optimizeFFT( dimension, dim, kernelSize, fft_planner_method, fftdimPerformed, lambdaCalib, flagZeropadding, sPrecision )
%OPTIMIZEFFT optimize the fft command
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

fprintf('Optimizing fft command\n');
if(strcmp(dimension,'2D'))
    if(dim(4) == 1)
        dimImg = [dim(1), dim(2), dim(5)]; % y-x-cha 
        dimKernel = [dimImg(1)+kernelSize(1)-1, dimImg(2)+kernelSize(2)-1, dimImg(3)];
    else
        dimImg = [dim(4), dim(1), dim(2), dim(5)]; % t-y-x-cha
        dimKernel = [dimImg(1), dimImg(2)+kernelSize(1)-1, dimImg(3)+kernelSize(2)-1, dimImg(4)];
    end    
elseif(strcmp(dimension,'3D'))
    dimImg = [dim(1), dim(3), dim(2), dim(5)]; % y-z-x-cha
    dimKernel = [dimImg(1)+kernelSize(1)-1, dimImg(2)+kernelSize(2)-1, dimImg(3)+kernelSize(3)-1, dimImg(4)];
elseif(strcmp(dimension,'4D'))
    dimImg = [dim(4), dim(1), dim(3), dim(2), dim(5)]; % t - y - z - x - cha
    dimKernel = [dimImg(1), dimImg(2)+kernelSize(1)-1, dimImg(3)+kernelSize(2)-1, dimImg(4)+kernelSize(3)-1, dimImg(5)];
elseif(strcmp(dimension,'5D'))
    dimImg = [dim(4), dim(1), dim(3), dim(2), dim(6), dim(5)]; % t - y - z - x - g - cha
    dimKernel = [dimImg(1), dimImg(2)+kernelSize(1)-1, dimImg(3)+kernelSize(2)-1, dimImg(4)+kernelSize(3)-1, dimImg(5), dimImg(6)];
end

% check fft along dimensions
% maximal fftdims = {1, 2, 3, [1 2], [1 3], [2 3], 1:3};
fftdims = {fftdimPerformed};
if(lambdaCalib > 0)
    fftdims{end+1} = 1:3; % due to corrKernelKspace -> convKernelImgspace
end

% check if wisdom file already exists
if(evalin('base','exist(''fftw_wisdom'',''var'')'))
    sizesCalculated = evalin('base', 'fftw_wisdom.sizesCalculated');
    doOptimization = false;
    if(size(sizesCalculated,3) ~= length(fftdims) || isemtpy(ismember(reshape(ipermute(sizesCalculated,[1 3 2]),2*size(sizesCalculated,3),size(sizesCalculated,2)),dimImg,'rows')) || isemtpy(ismember(reshape(ipermute(sizesCalculated,[1 3 2]),2*size(sizesCalculated,3),size(sizesCalculated,2)),dimKernel,'rows')))
        evalin('base', 'clear ''fftw_wisdom''');
        doOptimization = true;
    end
else
    doOptimization = true;
end
% if(exist([currpath,filesep,'utils',filesep,'general',filesep,'fftw_wisdom.mat'],'file'))
%     load([currpath,filesep,'utils',filesep,'general',filesep,'fftw_wisdom.mat']);
%     doOptimization = false;
%     if(size(sizesCalculated,3) ~= length(fftdims) || isemtpy(ismember(reshape(ipermute(sizesCalculated,[1 3 2]),2*size(sizesCalculated,3),size(sizesCalculated,2)),dimImg,'rows')) || isemtpy(ismember(reshape(ipermute(sizesCalculated,[1 3 2]),2*size(sizesCalculated,3),size(sizesCalculated,2)),dimKernel,'rows')))
%         delete([currpath,filesep,'utils',filesep,'general',filesep,'fftw_wisdom.mat']);
%         doOptimization = true;
%     end
% else
%     doOptimization = true;
% end

if(doOptimization)
    testImg = complex(randn(dimImg,sPrecision),randn(dimImg,sPrecision));
    if(flagZeropadding && lambdaCalib > 0)
        testKernel = complex(randn(dimKernel,sPrecision),randn(dimKernel,sPrecision));
    end
        
    sizesCalculated = zeros(2,length(dimImg),7);
    wisdom_str = cell(2,7);
    for i=1:length(fftdims)
        depth = sum(fftdims{i});
        if(length(fftdims{i}) > 1), depth = depth + 1; end;
        % optimize fft for image dimension
        fftw('planner',fft_planner_method);
        tmp = fftnshift(testImg,fftdims{i});
        wisdom_str{1,depth} = fftw('wisdom');
        sizesCalculated(1,:,depth) = dimImg;
        
        if(flagZeropadding && lambdaCalib > 0)
            % optimize fft for image kernel dimension
            fftw('planner',fft_planner_method);
            tmp = fftnshift(testKernel,fftdims{i});
            wisdom_str{2,depth} = fftw('wisdom');
            sizesCalculated(2,:,depth) = dimKernel;        
        end
    end  
    
    fftw_wisdom.sizesCalculated = sizesCalculated;
    fftw_wisdom.wisdom_str = wisdom_str;
    assignin('base', 'fftw_wisdom', fftw_wisdom);
%     save([currpath,filesep,'utils',filesep,'general',filesep,'fftw_wisdom.mat'], 'wisdom_str', 'sizesCalculated');
    clear 'testimg' 'testKernel';
end

