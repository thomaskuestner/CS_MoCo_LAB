function [fTrafo, bTrafo, kernelFTrafo, kernelBTrafo, para] = compBasis(kSpace,obj,dSensemap)
% compute transformation basis and provide for forward and backward
% transformation
% kSpace: 
% 2Dt: t  - k_y -  x  - cha  
% 2D: k_y - k_x - cha
% 3D: k_y - k_z -  x  - cha
% 4D:  t  - k_y - k_z -  x  - cha
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

if(nargin < 3)
    para.dSensemap = 1;
else
    para.dSensemap = dSensemap;
end
para.trafoType = obj.trafo.trafoType;
if(~isfield(obj.trafo,'fftdim'))
    para.fftdim = 1:2;
else
    para.fftdim = obj.trafo.fftdim;
end
if(~isfield(obj.trafo, 'scrambledim'))
    para.scrambledim = para.fftdim;
else
    para.scrambledim = obj.trafo.scrambledim;
end
% check input dimensionality
para.dimensionality = obj.measPara.dimension;
if(strcmp(para.dimensionality,'2D') && obj.measPara.dim(4) > 1)
    para.dimensionality = '2Dt';
end
para.zeroPad = obj.trafo.zeroPad;
para.rescaleInterp = obj.trafo.rescaleInterp;
para.size = size(kSpace);  % save size before transformation (and permutation)
if(strcmp(para.dimensionality,'4D'))
    if(length(para.size) < 5)
        tmp = para.size;
        para.size = ones(1,5);
        para.size(1:length(tmp)) = tmp;
    end
elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
    if(length(para.size) < 4)
        tmp = para.size;
        para.size = ones(1,4);
        para.size(1:length(tmp)) = tmp;
    end
elseif(strcmp(para.dimensionality,'2D'))  
    if(length(para.size) < 3)
        tmp = para.size;
        para.size = ones(1,3);
        para.size(1:length(tmp)) = tmp;
    end
end
para.shape = obj.trafo.shape; % 1D, 2D, 3D or 4D reconstruction
para.permRule = obj.trafo.permRule; % sparsifying dimensions are the first ones
if(~isempty(para.permRule)) % size before transformation, but after permutation
    para.imgsize = para.size(para.permRule); 
else
    para.imgsize = para.size;
end
para.kspaceTrafo = obj.trafo.kspaceTrafo; % kSpace transformation (for y and z)
para.precision = obj.measPara.precision;

switch para.trafoType
    case 'fft'
        para.windowing = obj.trafo.windowing;
        para.windowType = obj.trafo.windowType;
        para.windowOpt = obj.trafo.windowOpt;

    case 'pca'
        sigma = cell(1,para.shape);
        sigmaKernel = sigma;
        pc = sigma;
        pcKernel = sigma;
        if(obj.measPara.dim(5) == 1)
            padsize = zeros(1,ndims(kSpace));
        else
            padsize = zeros(1,ndims(kSpace)-1); % minus cha
        end

        if(length(padsize) == 2) % 2D
            helper = true(size(kSpace));
            m = [size(kSpace,2), size(kSpace,3), size(kSpace,4)];    
            % find the k-space center
            if(iscell(obj.calibSize))
                s = obj.calibSize{1}; % automatically determined calibration size
            else
                s = obj.calibSize; % predefined calibration size
            end
            idx = cell(1,length(s));
            for n=1:length(s)
                if(mod(s(n),2) == 0)
                    idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
                else
                    idx{n} = floor(m(n)/2)+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2)-1;
                end

                tmp = [idx{n}(1) <= 0, idx{n}(end) > m(n)];
                if(any(tmp))
                    if(all(tmp)), error('crop(): Both index out of bounds'); end;
                    hShift = [abs(idx{n}(1)) + 1, idx{n}(end) - m(n)];
                    op = {idx{n}(1), idx{n}(end); '+', '-'; '<', '>'; m(n) + 1, 0};
                    eval(sprintf('if(op{1,~tmp} %s hShift(tmp) %s %d), idx{n} = idx{n} %s hShift(tmp); else idx{n} = idx{n}(idx{n} %s %d);end;',op{2,tmp}, op{3,tmp}, op{4,tmp}, op{2,tmp}, op{3,~tmp}, op{4,~tmp}));
               end
            end
            helper(idx{1},idx{2},:) = false;
        elseif(length(padsize) == 3)
            helper = true(size(kSpace));
            m = size(kSpace);
            % find the k-space center
            s = obj.calibSize; 
            idx = cell(1,length(s));
            for n=1:length(s)
                if(mod(s(n),2) == 0)
                    idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
                else
                    idx{n} = floor(m(n)/2)+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2)-1;
                end

                tmp = [idx{n}(1) <= 0, idx{n}(end) > m(n)];
                if(any(tmp))
                    if(all(tmp)), error('crop(): Both index out of bounds'); end;
                    hShift = [abs(idx{n}(1)) + 1, idx{n}(end) - m(n)];
                    op = {idx{n}(1), idx{n}(end); '+', '-'; '<', '>'; m(n) + 1, 0};
                    eval(sprintf('if(op{1,~tmp} %s hShift(tmp) %s %d), idx{n} = idx{n} %s hShift(tmp); else idx{n} = idx{n}(idx{n} %s %d);end;',op{2,tmp}, op{3,tmp}, op{4,tmp}, op{2,tmp}, op{3,~tmp}, op{4,~tmp}));
               end
            end
            helper(idx{1},idx{2},idx{3},:) = false;
        elseif(length(padsize) == 4)
            helper = true(size(kSpace));
            m = [size(kSpace,2), size(kSpace,3), size(kSpace,4), size(kSpace,5)];
            if(iscell(obj.calibSize))
                % automatically determined calibration size
                loop = size(kSpace,1);
                pos = 't,idx{1},idx{2},idx{3},:';
            else
                % predefined calibration size
                loop = 1;
                pos = ':,idx{1},idx{2},idx{3},:';
            end
            for t=1:loop
                % find the k-space center
                if(iscell(obj.calibSize))
                    s = obj.calibSize{t};
                else
                    s = obj.calibSize;
                end

                idx = cell(1,length(s));
                for n=1:length(s)
                    if(mod(s(n),2) == 0)
                        idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
                    else
                        idx{n} = floor(m(n)/2)+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2)-1;
                    end

                    tmp = [idx{n}(1) <= 0, idx{n}(end) > m(n)];
                    if(any(tmp))
                        if(all(tmp)), error('crop(): Both index out of bounds'); end;
                        hShift = [abs(idx{n}(1)) + 1, idx{n}(end) - m(n)];
                        op = {idx{n}(1), idx{n}(end); '+', '-'; '<', '>'; m(n) + 1, 0};
                        eval(sprintf('if(op{1,~tmp} %s hShift(tmp) %s %d), idx{n} = idx{n} %s hShift(tmp); else idx{n} = idx{n}(idx{n} %s %d);end;',op{2,tmp}, op{3,tmp}, op{4,tmp}, op{2,tmp}, op{3,~tmp}, op{4,~tmp}));
                    end
                end
                eval(sprintf('helper(%s) = false;',pos));
            end 
        elseif(length(padsize) == 5)
            helper = true(size(kSpace));
            m = [size(kSpace,2), size(kSpace,3), size(kSpace,4), size(kSpace,5)];
            if(iscell(obj.calibSize))
                % automatically determined calibration size
                loop = size(kSpace,1);
                pos = 't,idx{1},idx{2},idx{3},:,:';
            else
                % predefined calibration size
                loop = 1;
                pos = ':,idx{1},idx{2},idx{3},:,:';
            end
            for t=1:loop
                % find the k-space center
                if(iscell(obj.calibSize))
                    s = obj.calibSize{t};
                else
                    s = obj.calibSize;
                end

                idx = cell(1,length(s));
                for n=1:length(s)
                    if(mod(s(n),2) == 0)
                        idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
                    else
                        idx{n} = floor(m(n)/2)+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2)-1;
                    end

                    tmp = [idx{n}(1) <= 0, idx{n}(end) > m(n)];
                    if(any(tmp))
                        if(all(tmp)), error('crop(): Both index out of bounds'); end;
                        hShift = [abs(idx{n}(1)) + 1, idx{n}(end) - m(n)];
                        op = {idx{n}(1), idx{n}(end); '+', '-'; '<', '>'; m(n) + 1, 0};
                        eval(sprintf('if(op{1,~tmp} %s hShift(tmp) %s %d), idx{n} = idx{n} %s hShift(tmp); else idx{n} = idx{n}(idx{n} %s %d);end;',op{2,tmp}, op{3,tmp}, op{4,tmp}, op{2,tmp}, op{3,~tmp}, op{4,~tmp}));
                    end
                end
                eval(sprintf('helper(%s) = false;',pos));
            end 
        end
        kSpace(helper) = 0; % just take low-res image

        if(~isempty(para.fftdim))
            tmp = ifftnshift(kSpace,para.fftdim,para.scrambledim); % y-z-x-cha
        end
        if(~isempty(para.permRule))
            tmp = permute(tmp,para.permRule);
        end 
        dim = size(tmp); % k_y-k_z-x-cha
        if(obj.measPara.dim(5) == 1)
            dim = [dim, 1];
        end
        if(obj.measPara.dim(4) > 1) % include t
            padsize(1) = para.size(1);
            if(obj.flagZeropadding)
                padsize(2:end) = para.size(2:end-1) + obj.kernelSize - 1;
            else
                padsize(2:end) = para.size(2:end-1);
            end
        else
            if(obj.flagZeropadding)
                padsize = para.size(1:end-1) + obj.kernelSize - 1;
            else
                padsize = para.size(1:end-1);
            end
        end
        if(~isempty(para.permRule))
            padsize = permute(padsize,para.permRule);
        end
%         para.permVec = [2 1 3 4 5; 1 2 3 4 5; 1 3 2 4 5]; % each row-n contains permutation vector so that n-th dimension is transformed
%         para.permVecTrafo = [1 2 3 4 5; 2 1 3 4 5; 3 1 2 4 5];

        if(length(padsize) == 2) % y-x(-cha)
            para.permVec = [2 1 3];
            para.permVecTrafo = [1 2 3];
        elseif(length(padsize) == 3) % t-y-x(-cha) OR y-z-x(-cha)
            para.permVec = [2 1 3 4; 1 2 3 4];
            para.permVecTrafo = [1 2 3 4; 2 1 3 4];
        elseif(length(padsize) == 4) % t-y-z-x(-cha)
            para.permVec = [2 3 1 4 5; 1 3 2 4 5; 1 2 3 4 5];
            para.permVecTrafo = para.permVec;
        else % t-y-z-g-x(-cha)
            para.permVec = [2 3 4 1 5 6; 1 3 4 2 5 6; 1 2 4 3 5 6];
            para.permVecTrafo = para.permVec;
        end

        for i=1:para.shape
%             sigma{i} = permute(tmp,para.permVec(i,:)); % TK

%             tmp = permute(tmp,para.permVec(i,:)); % TK
            % variante 1)
            if(length(padsize) == 2) % y-x(-cha)
               sigma{i} = permute(tmp,para.permVec(i,:));
            elseif(length(padsize) == 3) % t-y-x(-cha) OR y-z-x(-cha)
                sigma{i} = permute(tmp,para.permVec(i,:));
                sigma{i} = complex(zeros([size(tmp,para.permVec(i,1)), size(tmp,para.permVec(i,2)),  size(tmp,para.permVec(i,3)), size(tmp,para.permVec(i,4))],obj.measPara.precision),zeros([size(tmp,para.permVec(i,1)), size(tmp,para.permVec(i,2)),  size(tmp,para.permVec(i,3)), size(tmp,para.permVec(i,4))],obj.measPara.precision));
            elseif(length(padsize) == 4) % t-y-z-x(-cha)
                sigma{i} = complex(zeros([size(tmp,para.permVec(i,1))*size(tmp,para.permVec(i,2)),size(tmp,para.permVec(i,3)),size(tmp,para.permVec(i,4)),size(tmp,para.permVec(i,5))],obj.measPara.precision),zeros([size(tmp,para.permVec(i,1))*size(tmp,para.permVec(i,2)),size(tmp,para.permVec(i,3)),size(tmp,para.permVec(i,4)),size(tmp,para.permVec(i,5))],obj.measPara.precision));           
                for iCha=1:size(sigma{i},4)
                    for jDim=1:size(sigma{i},3)
    %                     sigma{i}(:,:,jDim,iCha) = reshape(permute(tmp(:,:,:,jDim,iCha),[1 3 2]),size(sigma{i},1),size(sigma{i},2));
                        sigma{i}(:,:,jDim,iCha) = reshape(tmp(:,:,:,jDim,iCha),size(sigma{i},1),size(sigma{i},2));
                    end
                end
            else % t-y-z-g-x(-cha)
                sigma{i} = complex(zeros([size(tmp,para.permVec(i,1))*size(tmp,para.permVec(i,2))*size(tmp,para.permVec(i,3)), size(tmp,para.permVec(i,4)), size(tmp,para.permVec(i,5)), size(tmp,para.permVec(i,6))],obj.measPara.precision),zeros([size(tmp,para.permVec(i,1))*size(tmp,para.permVec(i,2))*size(tmp,para.permVec(i,3)), size(tmp,para.permVec(i,4)), size(tmp,para.permVec(i,5)), size(tmp,para.permVec(i,6))],obj.measPara.precision));
                 for iCha=1:size(sigma{i},4)
                    for jDim=1:size(sigma{i},3)
    %                     sigma{i}(:,:,jDim,iCha) = reshape(permute(tmp(:,:,:,jDim,iCha),[1 3 2]),size(sigma{i},1),size(sigma{i},2));
                        sigma{i}(:,:,jDim,iCha) = reshape(tmp(:,:,:,:,jDim,iCha),size(sigma{i},1),size(sigma{i},2));
                    end
                end
            end


%             % variante 2)
%             sigma{i} = zeros([size(tmp,2)*size(tmp,3)*size(tmp,4),size(tmp,1),size(tmp,5)]); % tmp: t-y-z-x-cha
%             for iCha=1:size(sigma{i},4)
%                 sigma{i}(:,:,iCha) = reshape(permute(tmp(:,:,:,:,iCha),[1 3 2]),size(sigma{i},1),size(sigma{i},2));
%             end


            if(obj.flagZeropadding && obj.lambdaCalib > 0)
                if(para.zeroPad)
                    tmpKernel = zpad(tmp,[padsize(para.permVec(i,1:length(padsize))), nCha]);
                else
                    tmpKernel = complex(zeros([padsize(para.permVec(i,1:length(padsize))),nCha],obj.measPara.precision),zeros([padsize(para.permVec(i,1:length(padsize))),nCha],obj.measPara.precision));
                    if(length(padsize) == 2) % y-x(-cha)
                        for iCha=1:nCha
                            tmpKernel(:,:,iCha) = imresize(tmp(:,:,iCha), padsize(para.permVec(i,1:length(padsize))), para.rescaleInterp);
                        end
                    elseif(length(padsize) == 3) % t-y-x(-cha) OR y-z-x(-cha)
                        for iCha=1:nCha
                            RI = imref3d(size(tmp(:,:,:,iCha)),1,1,1);
                            scaleFactor = padsize(para.permVec(i,1:length(padsize)))./size(tmp(:,:,:,iCha));
                            tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                            tmpKernel(:,:,:,iCha) = crop(imwarp(tmp(:,:,:,iCha),RI,tFormResize,para.rescaleInterp),padsize(para.permVec(i,1:length(padsize))));
                        end
                    elseif(length(padsize) == 4) % t-y-z-x(-cha)
                        for iCha=1:nCha
                            for jDim=1:size(tmp,4)
                                RI = imref3d(size(tmp(:,:,:,jDim,iCha)),1,1,1);
                                scaleFactor = padsize(para.permVec(i,1:length(padsize)-1))./size(tmp(:,:,:,jDim,iCha));
                                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                                tmpKernel(:,:,:,jDim,iCha) = crop(imwarp(tmp(:,:,:,jDim,iCha),RI,tFormResize,para.rescaleInterp),padsizepara.permVec(i,1:length(padsize)-1));
                            end
                        end
                    else % t-y-z-g-x(-cha)
                        for iCha=1:nCha
                            for jDim=1:size(tmp,5)
                                for iSmp=1:size(tmp,4)
                                    RI = imref3d(size(tmp(:,:,:,iSmp,jDim,iCha)),1,1,1);
                                    scaleFactor = padsize(para.permVec(i,1:length(padsize)-1))./size(tmp(:,:,:,iSmp,jDim,iCha));
                                    tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                                    tmpKernel(:,:,:,iSmp,jDim,iCha) = crop(imwarp(tmp(:,:,:,iSmp,jDim,iCha),RI,tFormResize,para.rescaleInterp),padsizepara.permVec(i,1:length(padsize)-1));
                                end
                            end
                        end
                    end                  
                end
                if(length(padsize) == 2 || length(padsize) == 3) % y-x(-cha) || t-y-x(-cha) OR y-z-x(-cha)
                    sigmaKernel{i} = permute(tmpKernel,para.permVec(i,:));
                elseif(length(padsize) == 4) % t-y-z-x(-cha)
                    sigmaKernel{i} = complex(zeros([size(tmpKernel,para.permVec(i,1))*size(tmpKernel,para.permVec(i,2)),size(tmpKernel,para.permVec(i,3)),size(tmpKernel,para.permVec(i,4)),size(tmpKernel,para.permVec(i,5))],obj.measPara.precision),zeros([size(tmpKernel,para.permVec(i,1))*size(tmpKernel,para.permVec(i,2)),size(tmpKernel,para.permVec(i,3)),size(tmpKernel,para.permVec(i,4)),size(tmpKernel,para.permVec(i,5))],obj.measPara.precision));           
                    for iCha=1:size(sigmaKernel{i},4)
                        for jDim=1:size(sigmaKernel{i},3)
                              sigmaKernel{i}(:,:,jDim,iCha) = reshape(tmpKernel(:,:,:,jDim,iCha),size(sigmaKernel{i},1),size(sigmaKernel{i},2));
                        end
                    end
                else  % t-y-z-g-x(-cha)
                    sigmaKernel{i} = complex(zeros([size(tmpKernel,para.permVec(i,1))*size(tmpKernel,para.permVec(i,2))*size(tmpKernel,para.permVec(i,3)), size(tmpKernel,para.permVec(i,4)), size(tmpKernel,para.permVec(i,5)), size(tmpKernel,para.permVec(i,6))],obj.measPara.precision),zeros([size(tmpKernel,para.permVec(i,1))*size(tmpKernel,para.permVec(i,2))*size(tmpKernel,para.permVec(i,3)), size(tmpKernel,para.permVec(i,4)), size(tmpKernel,para.permVec(i,5)), size(tmpKernel,para.permVec(i,6))],obj.measPara.precision));
                    for iCha=1:size(sigmaKernel{i},4)
                        for jDim=1:size(sigmaKernel{i},3)
                            sigmaKernel{i}(:,:,jDim,iCha) = reshape(tmpKernel(:,:,:,:,jDim,iCha),size(sigmaKernel{i},1),size(sigmaKernel{i},2));
                        end
                    end                    
                end
%                 sigmaKernel{i} = permute(tmpKernel,para.permVec(i,:));
%                 sigmaKernel{i} = sigmaKernel{i} - repmat(mean(sigmaKernel{i}),[size(sigmaKernel{i},1),ones(1,length(padsize))]);
                pcKernel{i} = complex(zeros([padsize(para.permVec(i,1)),padsize(para.permVec(i,1)),padsize(para.permVec(i,2:length(padsize))),nCha],obj.measPara.precision),zeros([padsize(para.permVec(i,1)),padsize(para.permVec(i,1)),padsize(para.permVec(i,2:length(padsize))),nCha],obj.measPara.precision)); % nCha not included in padsize
            end
%             sigma{i} = sigma{i} - repmat(mean(sigma{i}),[size(sigma{i},1),ones(1,length(dim)-1)]);            
%             pc{i} = zeros([dim(para.permVec(i,2)),dim(para.permVec(i,2)),dim(para.permVec(i,3:length(dim)))]); % nCha included in dim
%             
            if(obj.flagZeropadding && obj.lambdaCalib > 0)
                itersize = [size(sigmaKernel{i},3),size(sigmaKernel{i},4)];
            else
                itersize = [size(sigma{i},3),size(sigma{i},4)];
            end
            if(length(padsize) == 2) % y-x(-cha)
                for iCha=1:dim(end)
%                     pc{i}(:,:,iCha) = pcacov(cov(sigma{i}(:,:,iCha))); % along para.permVec(i,1)-th dimension
                    try
                        pc{i}(:,:,iCha) = princomp(sigma{i}(:,:,iCha));
                    catch
                        pc{i}(:,:,iCha) = pca(sigma{i}(:,:,iCha));
                    end
                    if(obj.flagZeropadding && obj.lambdaCalib > 0)
                        try
                            pcKernel{i}(:,:,iCha) = princomp(sigmaKernel{i}(:,:,iCha)); % along para.permVec(i,1)-th dimension
                        catch
                            pcKernel{i}(:,:,iCha) = pca(sigmaKernel{i}(:,:,iCha));
                        end
                    end
                end
            else %if(length(padsize) == 3 || length(padsize) == 4 || length(padsize) == 5) % t-y-x(-cha) OR y-z-x(-cha) || t-y-z-x(-cha) || t-y-z-g-x(-cha)
                for iCha=1:dim(end)
                    for jDim=1:itersize(1)
                        if(jDim <= size(sigma{i},3))
                            try
                                pc{i}(:,:,jDim,iCha) = princomp(sigma{i}(:,:,jDim,iCha)); % along para.permVec(i,1)-th dimension
                            catch
                                pc{i}(:,:,jDim,iCha) = pca(sigma{i}(:,:,jDim,iCha)); % complete svd (slower)
                            end
                        end
                        if(obj.flagZeropadding && obj.lambdaCalib > 0)
                            try
                                pcKernel{i}(:,:,jDim,iCha) = princomp(sigmaKernel{i}(:,:,jDim,iCha)); % along para.permVec(i,1)-th dimension
                            catch
                                pcKernel{i}(:,:,jDim,iCha) = pca(sigmaKernel{i}(:,:,jDim,iCha));
                            end
                        end
                    end
                end
            end
        end
        para.pc = pc;
        para.pcKernel = pcKernel;
        para.permDim = dim;
        clear 'sigma' 'sigmaKernel' 'pc' 'pcKernel' 'tmp' 'tmpKernel';

    case 'wavelet_lab'
        para.waveletFilter = obj.trafo.waveletFilter;
        para.waveletFilterSize = obj.trafo.waveletFilterSize;    
        para.wavWeight = obj.trafo.wavWeight;
        if(~isempty(para.permRule))
            tmp = para.size(para.permRule);
        else
            tmp = para.size;
        end
        ssx = 2^ceil(log2(tmp(1))); 
        ssy = 2^ceil(log2(tmp(2)));
        para.ss = max(ssx, ssy);
        para.L = log2(para.ss/2^(obj.trafo.waveletStages)); % coarse level
        % create wavelet object
        para.W = Wavelet(obj.trafo.waveletFilter,obj.trafo.waveletFilterSize,para.L);
        para.flagThresholding = obj.trafo.flagThresholding;

    case 'wavelet_mat'
        para.waveletFilter = obj.trafo.waveletFilter;
        para.waveletFilterSize = obj.trafo.waveletFilterSize;
        if(strcmp(para.waveletFilter,'dmey'))
            para.wFilter = para.waveletFilter;
        else
            para.wFilter = [para.waveletFilter,num2str(para.waveletFilterSize)];
        end
        para.extMode = obj.trafo.extMode;
        para.waveletStages = obj.trafo.waveletStages;
        para.flagThresholding = obj.trafo.flagThresholding;
        % set fixed value for thresholding or estimate it from noise std
        if(isempty(obj.trafo.wavWeight))
            tmp = ifftnshift(kSpace,para.fftdim,para.scrambledim); % !!! or just central part
            tmp = tmp(:,:,:,1);
            noise = 1.4826 * median(abs(tmp(:) - median(tmp(:)))); %0.0056
            para.wavWeight = 3 * noise;
        else
            para.wavWeight = obj.trafo.wavWeight;
        end

    case 'mellin'
        para.sigma  	= obj.trafo.sigma; % sigma scaling value for AFMT
        para.extrapVal  = obj.trafo.extrapVal; % extrapolation value
        para.trafoSize  = obj.trafo.trafoSize; 
        para.interp     = obj.trafo.interp; % interpolation method
        para.bTrafoType = obj.trafo.bTrafoType; %type of backward trafo

    case 'curvelet'
        para.allCurvelets = obj.trafo.allCurvelets; % use wavelet or curvelet for finest scale (0/1)
        para.nbscales = obj.trafo.nbscales; %floor(log2(min([para.size(1),para.size(2),para.size(3)])))-2;  % number of radial scales (calculated as given in the example fdct3d_demo_basic.m)
        para.nbdstz_coarse = obj.trafo.nbdstz_coarse; % number of angular scales (has to be a multiple of 4)

    case 'surfacelet'
        para.Pyr_mode = obj.trafo.Pyr_mode;
        para.HGfname = obj.trafo.HGfname;
        para.bo = obj.trafo.bo;
        para.decompLevels = obj.trafo.decompLevels;
        para.msize = obj.trafo.msize;
        para.beta = obj.trafo.beta;
        para.lambda = obj.trafo.lambda;
        para.Lev_array = obj.trafo.Lev_array;
        para.downSamp = obj.trafo.downSamp;       
        if(obj.trafo.padOnSameSize && obj.flagZeropadding)
            maxsize = para.size(1:length(obj.kernelSize)) + obj.kernelSize - 1; % y-z-x OR y-x
            if(~isempty(para.permRule))
                maxsize = maxsize(para.permRule);
            end
            maxsize = max(maxsize(1:para.shape)); % maximal size occurs for kernelImg
            n = maxsize./para.downSamp;
            n(n<1) = 1;
            n = ceil(n);
            para.padsize = max(n.*para.downSamp);
        else
            para.padsize = [];
        end

%         para.padSize = lcm(maxsize,2^para.decompLevels);
%         if(para.Pyr_mode == 1)
%             para.padSize = max(para.size);
%         else
%             tmp = para.Pyr_mode^(para.decompLevels - 1) * para.Pyr_mode;
%             tmp = tmp:tmp:max(para.size)+tmp;
%             para.padSize = tmp(end);
%         end  

    case 'bandlet'

    case 'gabor'
        % http://stackoverflow.com/questions/9003147/how-to-apply-gabor-wavelets-to-an-image

    case 'wavelet_packet' % complex dual-tree 3D wavelet

end

% free memory
clear 'kSpace'

% positions
% 1: initialization of W    cart => new_basis & IFFT2 = forwardTransform
% 2: error e                new_basis => cart & FFT2  = backwardTransform
% 3: gradient G             cart => new_basis & IFFT2 = forwardTransform
% 4: Z matrix               new_basis => cart & FFT2  = backwardTransform
% 5: dImg                   new_basis => cart         = outTransform

fTrafo = @forwardTransform;  % k_y-k_z-x-cha (cart)  => y-z-x-cha (new_basis) (IFFT-in wrapper)
bTrafo = @backwardTransform; % y-z-x-cha (new_basis) => k_y-k_z-x-cha (cart)  (FFT-out wrapper)
kernelFTrafo = @kernelFTransform; % y-z-x-cha (cart)      => y-z-x-cha (new_basis)
kernelBTrafo = @kernelBTransform; % y-z-x-cha (new_basis) => y-z-x-cha (cart)

if(obj.trafo.kspaceTrafo)
    % directly transform kSpace data
    fTrafo = @kernelFTransform; % k_y-k_z-x-cha (cart)      => k_y-k_z-x-cha (new_basis)
    bTrafo = @kernelBTransform; % k_y-k_z-x-cha (new_basis) => k_y-k_z-x-cha (cart)
end


    function img = forwardTransform(img)
        % k_y-k_z-x-cha (cart) => y-z-x-cha (new_basis)
        img = kernelFTransform(ifftnshift(img,para.fftdim,para.scrambledim).*conj(para.dSensemap));             
    end


    function img = kernelFTransform(img)   
        % permute image for sparsifying dimensions to be the preceding ones
        if(para.kspaceTrafo && ~isempty(para.fftdim)) % still transform remaining dimensions which are not k-space transformed
            img = ifftnshift(img,para.fftdim,para.scrambledim);
        end
        if(~isempty(para.permRule))
            img = permute(img,para.permRule);
        end       
        
        % y-z-x-cha (cart) => y-z-x-cha (new_basis)
        switch para.trafoType
            case 'fft'
                if(para.windowing)
                    if(para.shape == 1)
                        [~,window] = windowND(para.windowType,size(img,1),para.windowOpt{1},para.windowOpt{2});
                        window = repmat(window,[1 size(img,2) size(img,3) size(img,4)]);
                    elseif(para.shape == 2)
                        window = windowND(para.windowType,[size(img,1) size(img,2)],para.windowOpt{1},para.windowOpt{2});
                        window = repmat(window,[1 1 size(img,3) size(img,4)]);
                    else
                        window = windowND(para.windowType,[size(img,1) size(img,2) size(img,3)],para.windowOpt{1},para.windowOpt{2});
                        window = repmat(window,[1 1 1 size(img,4)]);
                    end
                    img = img .* window;
                end
                clear 'window';
                
            case 'dct'
                % cart => DCT
                if(para.shape == 1)
                    for i=1:size(img,5)
                        for j=1:size(img,4)
                            for k=1:size(img,3)
                                for l=1:size(img,2)
                                    img(:,l,k,j,i) = dct(img(:,l,k,j,i));
                                end
                            end
                        end
                    end                    
                elseif(para.shape == 2)
                    for i=1:size(img,5)
                        for j=1:size(img,4)
                            for k=1:size(img,3)
                                img(:,:,k,j,i) = dct2(img(:,:,k,j,i));
                            end
                        end
                    end
                else
                    for lDCT = 1:length(para.shape)
                        permRule = 1:ndims(img);
                        permRule(2:end) = permRule(~ismember(permRule,para.shape(lDCT)));
                        permRule(1) = para.shape(lDCT);
                        img = permute(img,permRule);

                        for i=1:size(img,2)
                            for j=1:size(img,3)
                                for k=1:size(img,4)
                                    for l=1:size(img,5)
                                        img(:,i,j,k,l) = dct(img(:,i,j,k,l));
                                    end
                                end
                            end
                        end
                        img = ipermute(img,permRule);
                    end
                end
                                
            case 'pca'
                % cart => PCA
%                 n_dim = ndims(img);
                imgsize = size(img);
                for i=1:para.shape
                    img = permute(img,para.permVecTrafo(i,:));
                    imgsizePerm = size(img);
                    if(length(imgsize) == 5)
%                         img = permute(reshape(img,[imgsize(2)*imgsize(3), imgsize(1), imgsize(4), imgsize(5)]),[2 1 3 4]);
                        img = permute(reshape(img,[size(img,1)*size(img,2), size(img,3), size(img,4), size(img,5)]),[2 1 3 4]);
                    elseif(length(imgsize) == 6)
                        img = permute(reshape(img,[size(img,1)*size(img,2)*size(img,3), size(img,4), size(img,5), size(img,6)]),[2 1 3 4]);
                    end
                    img = mtimesx(para.pc{i},'C',img);
                    if(length(imgsize) == 5)
                        img = reshape(permute(img,[2 1 3 4]),imgsizePerm(1),imgsizePerm(2),imgsizePerm(3),imgsizePerm(4),imgsizePerm(5));
                    elseif(length(imgsize) == 6)
                        img = reshape(permute(img,[2 1 3 4]),imgsizePerm(1),imgsizePerm(2),imgsizePerm(3),imgsizePerm(4),imgsizePerm(5),imgsizePerm(6));
                    end
                    img = ipermute(img, para.permVecTrafo(i,:));
                    
%                     % t-y-z-x-cha
%                     img = permute(img,[2 3 1 4 5]); % y-z-t-x-cha
%                     tmp = permute(reshape(img,[imgsize(2)*imgsize(3), imgsize(1), imgsize(4), imgsize(5)]),[2 1 3 4]);
%                     tmp = permute(mtimesx(para.pc{i},'C',tmp),[2 1 3 4]); % t - y*z - x - cha
%                     img = permute(reshape(tmp,[imgsize(2), imgsize(3), imgsize(1), imgsize(4), imgsize(5)]),[3 1 2 4 5]);
                end
                
            case 'wavelet_lab'
                % cart => wavelet
                % just 2D shape, but arbitrary input dimensionality (2D-4D)
                imgsize = ones(1,5);
                imgTmp = size(img);
                imgsize(1:length(imgTmp)) = imgTmp;
                if(para.zeroPad)
                    img = zpad(img,[para.ss,para.ss,imgsize(3:end)]);
                else
                    imgOut = zeros(para.ss,para.ss,imgsize(3:end),para.precision);
                    for i=1:imgsize(5)
                        for j=1:imgsize(4)
                            for k=1:imgsize(3)
                                imgOut(:,:,k,j,i) = imresize(img(:,:,k,j,i), [para.ss,para.ss], para.rescaleInterp);
                            end
                        end
                    end
                    img = imgOut;
                    clear 'imgOut';
                end
                for i=1:imgsize(5)
                    for j=1:imgsize(4)
                        for k=1:imgsize(3)
                            img(:,:,k,j,i) = para.W * img(:,:,k,j,i);
                        end
                    end
                end
                % soft-thresholding
                if(para.flagThresholding)
                    img = SoftThresh(img,para.wavWeight);
                end
                
            case 'wavelet_mat'
                % cart => wavelet
                imgsize = size(img);
                imgIn = img;
                img = cell(1,para.imgsize(end));
                recinfo = img;
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:size(imgIn,5) % cha
                        if(para.shape == 3)
                            img{h} = cell(7*para.waveletStages+1,size(imgIn,4));
                            recinfo{h} = cell(1,size(imgIn,4));
                        end
                        for i=1:size(imgIn,4)
                            if(para.shape == 3)
                                helper = wavedec3(imgIn(:,:,:,i,h),para.waveletStages,para.wFilter,'mode',para.extMode);
                                img{h}(:,i) = helper.dec;
                                helper = rmfield(helper,'dec');
                                recinfo{h} = helper;
                                dispProgress('Trafo', ((h-1)*size(imgIn,4)*size(imgIn,3)*size(imgIn,2) + i*size(imgIn,3)*size(imgIn,2))/prod(imgsize(2:end)));
                                continue;
                            elseif(para.shape == 2)
                                img{h} = zeros(para.waveletStages+2,2,size(imgIn,3),size(imgIn,4),para.precision);
                                recinfo{h} = img{h};
                            end
                            for j=1:size(imgIn,3) % freq
                                if(para.shape == 2)
                                    [img{h}(:,:,j,i), recinfo{h}(:,:,j,i)] = wavedec2(imgIn(:,:,j,i,h),para.waveletStages,para.wFilter,'mode',para.extMode);
                                    dispProgress('Trafo', ((h-1)*size(imgIn,4)*size(imgIn,3)*size(imgIn,2) + (i-1)*size(imgIn,3)*size(imgIn,2) + j*size(imgIn,2))/prod(imgsize(2:end)));
                                    continue;
                                elseif(para.shape == 1)
                                    img{h} = zeros(para.waveletStages+2,size(imgIn,2),size(imgIn,3),size(imgIn,4),para.precision);
                                    recinfo{h} = img{h};
                                end
                                for k=1:size(imgIn,2) % z
                                    if(para.shape == 1)
                                        [img{h}(:,k,j,i),  recinfo{h}(:,k,j,i)] = wavedec(imgIn(:,k,j,i,h),para.waveletStages,para.wFilter,'mode',para.extMode);
                                    end
                                    dispProgress('Trafo', ((h-1)*size(imgIn,4)*size(imgIn,3)*size(imgIn,2) + (i-1)*size(imgIn,3)*size(imgIn,2) + (j-1)*size(imgIn,2) + k)/prod(imgsize(2:end)));
                                end
                                dispProgress('Trafo', ((h-1)*size(imgIn,4)*size(imgIn,3)*size(imgIn,2) + (i-1)*size(imgIn,3)*size(imgIn,2) + j*size(imgIn,2))/prod(imgsize(2:end)));
                            end
                            clear 'helper';
                            dispProgress('Trafo', ((h-1)*size(imgIn,4)*size(imgIn,3)*size(imgIn,2) + i*size(imgIn,3)*size(imgIn,2))/prod(imgsize(2:end)));
                        end
                        dispProgress('Trafo', (h*size(imgIn,4)*size(imgIn,3)*size(imgIn,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:size(imgIn,4) % cha
                        if(para.shape == 3)
                            helper = wavedec3(imgIn(:,:,:,i),para.waveletStages,para.wFilter,'mode',para.extMode);
                            img{i} = helper.dec;
                            helper = rmfield(helper,'dec');
                            recinfo{i} = helper;
                            dispProgress('Trafo', (i*size(imgIn,3)*size(imgIn,2))/prod(imgsize(2:end)));
                            continue;
                        elseif(para.shape == 2)
%                             img{i} = zeros(1,3*para.waveletStages+1,size(imgIn,3));
                            recinfo{i} = zeros(para.waveletStages+2,2,size(imgIn,3),para.precision);
                        end
                        for j=1:size(imgIn,3) % freq
                            if(para.shape == 2)
                                [helper, recinfo{i}(:,:,j)] = wavedec2(imgIn(:,:,j,i),para.waveletStages,para.wFilter);
                                if(j==1)
                                    img{i} = zeros([size(helper),size(imgIn,3)],para.precision);
                                end
                                img{i}(:,:,j) = helper;
                                dispProgress('Trafo', ((i-1)*size(imgIn,3)*size(imgIn,2) + j*size(imgIn,2))/prod(imgsize(2:end)));
                                continue;
                            elseif(para.shape == 1)
                                img{i} = zeros(para.waveletStages+2,size(imgIn,2),size(imgIn,3),para.precision);
                                recinfo{i} = img{i};
                            end
                            for k=1:size(imgIn,2) % z
                                if(para.shape == 1)
                                    [img{i}(:,k,j),  recinfo{i}(:,k,j)] = wavedec(imgIn(:,k,j,i),para.waveletStages,para.wFilter,'mode',para.extMode);
                                end
                                dispProgress('Trafo', ((i-1)*size(imgIn,3)*size(imgIn,2) + (j-1)*size(imgIn,2) + k)/prod(imgsize(2:end)));
                            end
                            dispProgress('Trafo', ((i-1)*size(imgIn,3)*size(imgIn,2) + j*size(imgIn,2))/prod(imgsize(2:end)));
                        end
                        clear 'helper';
                        dispProgress('Trafo', (i*size(imgIn,3)*size(imgIn,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:size(imgIn,3) % cha
                        if(para.shape == 2)
                            [img{i}, recinfo{i}] = wavedec2(imgIn(:,:,i),para.waveletStages,para.wFilter,'mode',para.extMode);
                            dispProgress('Trafo', (i*size(imgIn,2))/prod(imgsize(2:end)));
                            continue;
                        elseif(para.shape == 1)
                            img{i} = zeros(para.waveletStages+2,size(imgIn,2),para.precision);
                            recinfo{i} = img{i};
                        end
                        for j=1:size(imgIn,2)
                            if(para.shape == 1)
                                [img{i}(:,j), recinfo{i}(:,j)] = wavedec(imgIn(:,j,i),para.waveletStages,para.wFilter);
                                dispProgress('Trafo', ((i-1)*size(imgIn,2) + j)/prod(imgsize(2:end)));
                            end
                        end
                        dispProgress('Trafo', (i*size(imgIn,2))/prod(imgsize(2:end)));
                    end
                end
                dispProgress('Trafo', 'Close');
                % soft-thresholding
                if(para.flagThresholding)
                    for iCha=1:length(img)
                        if(iscell(img{iCha}))
                            for iJ=1:length(img{iCha})
                                img{iCha}{iJ} = SoftThresh(img{iCha}{iJ},para.wavWeight);
                            end
                        else          
                            img{iCha} = SoftThresh(img{iCha},para.wavWeight);
                        end
                    end
%                     absImg = img{1};
%                     for i=1:numel(img{1})
%                         tmp = cellfun(@(x) x{i}, img(2:nCha), 'UniformOutput', false);
%                         absImg{i}(:,:,:,2:nCha) = cell2mat(shiftdim(tmp,-2));
%                     end
%                     absImg = cellfun(@(x) sqrt(sum(abs(x).^2,4)), absImg, 'UniformOutput', false);
%                     res = cellfun(@(x) x-para.wavWeight, absImg, 'UniformOutput', false);
%                     res = cellfun(@(x) (x + abs(x))/2, res, 'UniformOutput', false);
%                     for iCha=1:nCha
%                         unity = cellfun(@(x,y) x./(y+eps), img{iCha}, absImg, 'UniformOutput', false);
%                         img{iCha} = cellfun(@(x,y) x.*y, unity, res, 'UniformOutput', false);
%                     end
                end
                img = TRAFO(img,recinfo);
                
            case 'mellin'
                % cart => mellin
                % in: y - z - x - cha (cart)
                % out: y - z - x - cha (mellin)
                imgsize = size(img);
                if(para.shape == 2)
                    dim = [size(img,2)*para.trafoSize(2), size(img,1)*para.trafoSize(1)];%2D trafo allows arbitrary transformation size
                    imgTra = cell(1, imgsize(end));
                    [imgTra{:}] = deal(zeros(dim(2), dim(1), imgsize(3:end-1),para.precision));
                end
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:size(img,5) % cha
                        for i=1:size(img,4)
                            if(para.shape == 3)
                                img(:,:,:,i,h) = AFMT3(img(:,:,:,i,h), para.sigma, size(img,2), size(img,1), size(img,3), para.extrapVal, para.interp, 'F-AFMT');
                                dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                                continue;
                            end
                            for j=1:size(img,3) % freq
                                if(para.shape == 2)
                                    imgTra{h}(:,:,j,i) = AFMT2(img(:,:,j,i,h), para.sigma, dim(1), dim(2), para.extrapVal, para.interp, 'F-AFMT');
                                    dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                                end
                            end
                            dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                        end
                        dispProgress('Trafo', (h*size(img,4)*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:size(img,4) % cha
                        if(para.shape == 3)
                            img(:,:,:,i) = AFMT3(img(:,:,:,i), para.sigma, size(img,2), size(img,1), size(img,3), para.extrapVal, para.interp, 'F-AFMT');
                            dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                            continue;
                        end
                        for j=1:size(img,3) % freq
                            if(para.shape == 2)
                                imgTra{i}(:,:,j) = AFMT2(img(:,:,j,i), para.sigma, dim(1), dim(2), para.extrapVal, para.interp, 'F-AFMT');
                                dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                            end
                        end
                        dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:size(img,3)
                        if(para.shape == 2)
                            imgTra{i}(:,:) = AFMT2(img(:,:,i), para.sigma, dim(1), dim(2), para.extrapVal, para.interp, 'F-AFMT');
                        end
                        dispProgress('Trafo', (i*size(img,2))/prod(imgsize(2:end)));
                    end
                end
                dispProgress('Trafo', 'Close');
                if(para.shape == 2)
                    img = TRAFO(imgTra, size(img)); 
                    clear 'imgTra';
                end
                
            case 'curvelet'
                % cart => curvelet
                imgsize = size(img);                
                imgTra = cell(1, imgsize(end));
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:size(img,5) % cha
                        if(para.shape == 3)
                            imgTra{h} = cell(size(img,4),1);
                        elseif(para.shape == 2)
                            imgTra{h} = cell(size(img,4),size(img,3));
                        end
                        for i=1:size(img,4)
                            if(para.shape == 3)
                                imgTra{h}{i} =  fdct3d_forward_mex(size(img, 1), size(img, 2), size(img, 3), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,:,i,h));                                
                                dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                                continue;
                            end
                            for j=1:size(img,3) % freq
                                if(para.shape == 2)
                                    imgTra{h}{i,j} = fdct_usfft_mex(size(img,1), size(img,2), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,j,i,h));
                                    dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                                end
                            end
                            dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                        end
                        dispProgress('Trafo', (h*size(img,4)*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:size(img,4) % cha
                        if(para.shape == 3)
                            imgTra{i} = fdct3d_forward_mex(size(img, 1), size(img, 2), size(img, 3), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,:,i));
                            dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                            continue;
                        elseif(para.shape == 2)
                            imgTra{i} = cell(size(img,3),1);
                            for j=1:size(img,3)
                               imgTra{i}{j} = fdct_usfft_mex(size(img,1), size(img,2), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,j,i));
                               dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                            end
                        end
                        dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:size(img,3)
                        if(para.shape == 2)
                            imgTra{i} = fdct_usfft_mex(size(img,1), size(img,2), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,i));
                            dispProgress('Trafo', (i*size(img,2))/prod(imgsize(2:end)));
                        end     
                    end
                end
                dispProgress('Trafo', 'Close');
                img = TRAFO(imgTra, size(img)); %store original size of img as meta data
                    
            case 'surfacelet'
                % cart => surfacelet
                nCha = para.size(end);
                imgsize = size(img);
                if(isempty(para.padsize))
                    padsize = imgsize; % already permuted image
    %                 n = ceil(padsize(1:3)/2^para.decompLevels);
    %                 padsize = max(n.*2^para.decompLevels);
                    n = padsize(1:para.shape)./para.downSamp;
                    n(n<1) = 1;
                    n = ceil(n);
                    padsize = max(n.*para.downSamp);
                else
                    padsize = para.padsize;
                end
                para.padsize = padsize; % size of imgIn
                % para.imgsize = imgsize; % size of img (after permutation)
                if(para.zeroPad)
                    imgIn = zpad(img,[padsize*ones(1,para.shape),imgsize(para.shape+1:end)]);
                else
                    imgIn = zeros([padsize*ones(1,para.shape),imgsize(para.shape+1:end)],para.precision);
                    if(para.shape == 3) % 3D padding
                        for i=1:size(img,5)
                            for j=1:size(img,4)
                                RI = imref3d(size(img(:,:,:,j,i)),1,1,1);
                                scaleFactor = size(imgIn(:,:,:,j,i))./size(img(:,:,:,j,i));
                                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                                imgIn(:,:,:,j,i) = crop(imwarp(img(:,:,:,j,i),RI,tFormResize,para.rescaleInterp),padsize*ones(1,para.shape));
                            end
                        end                        
                    elseif(para.shape == 2)
                        for i=1:size(img,5)
                            for j=1:size(img,4)
                                for k=1:size(img,3)
                                    imgIn(:,:,k,j,i) = imresize(img(:,:,k,j,i), [padsize padsize], para.rescaleInterp);
                                end
                            end
                        end
                    end
                end   
                clear 'img'
                img = cell(1,nCha);
%                 recinfo = cell(1,nCha);
                recinfo = cell(1,2); % reduce overhead
                imgsize = size(imgIn);
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:size(imgIn,5) % cha
                        if(para.shape == 3)
                            img{h} = cell(1,size(imgIn,4));
%                             recinfo{h} = img{h};
                        end
                        for i=1:size(imgIn,4)
                            if(para.shape == 3)
                                [imgReal, recinfo{1,1}] = Surfdec(real(imgIn(:,:,:,i,h)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda); % recinfo{1,h}{1,i}
                                [imgImag, ~] = Surfdec(imag(imgIn(:,:,:,i,h)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda);
                                [out, idx] = flattenCellMatrix(imgReal);
                                recinfo{1,2} = idx;
                                clear 'imgReal'
                                [outB, idxB] = flattenCellMatrix(imgImag);
                                clear 'imgImag'
                                if(~all(cellfun(@(x,y) isequal(x,y), idx, idxB)))
                                    error('compBasis: Nested cells must have the same size');
                                end
                                img{1,h}{1,i} = cellfun(@(x,y) x+1i*y, out, outB, 'UniformOutput', false);
                                clear 'out' 'outB' 'idxB'
%                                 img{1,h}{1,i} = reconFlatCellMatrix(img{1,h}{1,i},idx);
                                dispProgress('Trafo', ((h-1)*size(imgIn,4)*size(imgIn,3)*size(imgIn,2) + i*size(imgIn,3)*size(imgIn,2))/prod(imgsize(2:end)));
                                continue;
                            elseif(para.shape == 2)
%                                 img{1,h}{1,i} = cell(length(para.Lev_array)+1,size(imgIn,3));
%                                 recinfo{1,h}{1,i} = img{1,h}{1,i};
                            end
                            for j=1:size(imgIn,3)
                                if(para.shape == 2)
                                    [imgReal, recinfo{1,1}] = Surfdec(real(imgIn(:,:,j,i,h)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda); % recinfo{1,h}{1,i}
                                    [imgImag, ~] = Surfdec(imag(imgIn(:,:,:,j,i,h)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda);
                                    [out, idx] = flattenCellMatrix(imgReal);
                                    recinfo{1,2} = idx;
                                    clear 'imgReal'
                                    [outB, idxB] = flattenCellMatrix(imgImag);
                                    clear 'imgImag'
                                    if(~all(cellfun(@(x,y) isequal(x,y), idx, idxB)))
                                        error('compBasis: Nested cells must have the same size');
                                    end
                                    if(j==1), img{1,h}{1,i} = cell([size(out), size(imgIn,3)]); end;
                                    img{1,h}{1,i}(:,:,j) = cellfun(@(x,y) x+1i*y, out, outB, 'UniformOutput', false);
                                    clear 'out' 'outB' 'idxB'
%                                     img{1,h}{1,i}(:,j) = reconFlatCellMatrix(tmp,idx);
%                                     clear 'tmp'
                                    dispProgress('Trafo', ((h-1)*size(imgIn,4)*size(imgIn,3)*size(imgIn,2) + (i-1)*size(imgIn,3)*size(imgIn,2) + j*size(imgIn,2))/prod(imgsize(2:end)));
                                end
                            end
                        end
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:size(imgIn,4) % cha
                        if(para.shape == 3)
                            [imgReal, recinfo{1,1}] = Surfdec(real(imgIn(:,:,:,i)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda); % recinfo{1,i}
                            [imgImag, ~] = Surfdec(imag(imgIn(:,:,:,i)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda);
                            [out, idx] = flattenCellMatrix(imgReal);
                            recinfo{1,2} = idx;
                            clear 'imgReal'
                            [outB, idxB] = flattenCellMatrix(imgImag);
                            clear 'imgImag'
                            if(~all(cellfun(@(x,y) isequal(x,y), idx, idxB)))
                                error('compBasis: Nested cells must have the same size');
                            end
                            img{1,i} = cellfun(@(x,y) x+1i*y, out, outB, 'UniformOutput', false);
                            clear 'out' 'outB' 'idxB'
%                             img{1,i} = reconFlatCellMatrix(img{1,i},idx);
                            dispProgress('Trafo', (i*size(imgIn,3)*size(imgIn,2))/prod(imgsize(2:end)));
                            continue;
                        elseif(para.shape == 2)
%                             img{1,i} = cell(length(para.Lev_array)+1,size(imgIn,3));
%                             recinfo{1,i} = img{1,i};
                        end
                        for j=1:size(imgIn,3)
                            if(para.shape == 2)
                                [imgReal, recinfo{1,1}] = Surfdec(real(imgIn(:,:,j,i)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda); % recinfo{1,i}
                                [imgImag, ~] = Surfdec(imag(imgIn(:,:,j,i)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda);
                                [out, idx] = flattenCellMatrix(imgReal);
                                recinfo{1,2} = idx;
                                clear 'imgReal'
                                [outB, idxB] = flattenCellMatrix(imgImag);
                                clear 'imgImag'
                                if(~all(cellfun(@(x,y) isequal(x,y), idx, idxB)))
                                    error('compBasis: Nested cells must have the same size');
                                end
                                if(j==1), img{1,i} = cell([size(out), size(imgIn,3)]); end;
                                img{1,i}(:,:,j) = cellfun(@(x,y) x+1i*y, out, outB, 'UniformOutput', false);
                                clear 'out' 'outB' 'idxB'
%                                 img{1,i}(:,j) = reconFlatCellMatrix(tmp,idx);
%                                 clear 'tmp'
                                dispProgress('Trafo', ((i-1)*size(imgIn,3)*size(imgIn,2) + j*size(imgIn,2))/prod(imgsize(2:end)));
                            end
                        end
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:size(imgIn,3) % cha
                        if(para.shape == 2)
                            [imgReal, recinfo{1,1}] = Surfdec(real(imgIn(:,:,i)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda); % recinfo{1,i}
                            [imgImag, ~] = Surfdec(imag(imgIn(:,:,i)), para.Pyr_mode, para.Lev_array, para.HGfname, 'bo', para.bo, 'msize', para.msize, 'beta', para.beta, 'lambda', para.lambda);
                            [out, idx] = flattenCellMatrix(imgReal);
                            recinfo{1,2} = idx;
                            clear 'imgReal'
                            [outB, idxB] = flattenCellMatrix(imgImag);
                            clear 'imgImag'
                            if(~all(cellfun(@(x,y) isequal(x,y), idx, idxB)))
                                error('compBasis: Nested cells must have the same size');
                            end
                            img{1,i} = cellfun(@(x,y) x+1i*y, out, outB, 'UniformOutput', false);
                            clear 'out' 'outB' 'idxB'
%                             img{1,i} = reconFlatCellMatrix(img{1,i},idx);
                            dispProgress('Trafo', (i*size(imgIn,2))/prod(imgsize(2:end)));
                        end
                    end
                end
                dispProgress('Trafo', 'Close');
                clear 'imgIn'
                meta.recinfo = recinfo;
                img = TRAFO(img,meta);
        end        
    end


    function img = backwardTransform(img)
        % y-z-x-cha (new_basis) => k_y-k_z-x-cha (cart)
        img = fftnshift(kernelBTransform(img).*para.dSensemap,para.fftdim,para.scrambledim);
    end


    function img = kernelBTransform(img)
        % y-z-x-cha (new_basis) => y-z-x-cha (cart)
        switch para.trafoType
            case 'fft'
                if(para.windowing)
                    if(para.shape == 1)
                        [~,window] = windowND(para.windowType,size(img,1),para.windowOpt{1},para.windowOpt{2});
                        window = repmat(window,[1 size(img,2) size(img,3) size(img,4)]);
                    elseif(para.shape == 2)
                        window = windowND(para.windowType,[size(img,1) size(img,2)],para.windowOpt{1},para.windowOpt{2});
                        window = repmat(window,[1 1 size(img,3) size(img,4)]);
                    else
                        window = windowND(para.windowType,[size(img,1) size(img,2) size(img,3)],para.windowOpt{1},para.windowOpt{2});
                        window = repmat(window,[1 1 1 size(img,4)]);
                    end
                    img = img .* window;
                end
                clear 'window';

            case 'dct'
                % DCT => cart
%                 img = double(img);
                if(para.shape == 1)
                    for i=1:size(img,5)
                        for j=1:size(img,4)
                            for k=1:size(img,3)
                                for l=1:size(img,2)
                                    img(:,l,k,j,i) = idct(img(:,l,k,j,i));
                                end
                            end
                        end
                    end
                elseif(para.shape == 2)
                    for i=1:size(img,5)
                        for j=1:size(img,4)
                            for k=1:size(img,3)
                                img(:,:,k,j,i) = idct2(img(:,:,k,j,i));
                            end
                        end
                    end
                else
                    for lDCT = 1:length(para.shape)
                        permRule = 1:ndims(img);
                        permRule(2:end) = permRule(~ismember(permRule,para.shape(lDCT)));
                        permRule(1) = para.shape(lDCT);
                        img = permute(img,permRule);

                        for i=1:size(img,2)
                            for j=1:size(img,3)
                                for k=1:size(img,4)
                                    for l=1:size(img,5)
                                        img(:,i,j,k,l) = idct(img(:,i,j,k,l));
                                    end
                                end
                            end
                        end
                        img = ipermute(img,permRule);
                    end
                end
                
            case 'pca'
                % PCA => cart
%                 n_dim = ndims(img);
                imgsize = size(img);
                for i=1:para.shape
                    img = permute(img,para.permVecTrafo(i,:));
                    imgsizePerm = size(img);
                    if(length(imgsize) == 5)
%                         img = permute(reshape(img,[imgsize(2)*imgsize(3), imgsize(1), imgsize(4), imgsize(5)]),[2 1 3 4]);
                        img = permute(reshape(img,[size(img,1)*size(img,2), size(img,3), size(img,4), size(img,5)]),[2 1 3 4]);
                    elseif(length(imgsize) == 6)
                        img = permute(reshape(img,[size(img,1)*size(img,2)*size(img,3), size(img,4), size(img,5), size(img,6)]),[2 1 3 4]);
                    end
                    img = mtimesx(para.pc{i},img);
                    if(length(imgsize) == 5)
                        img = reshape(permute(img,[2 1 3 4]),imgsizePerm(1),imgsizePerm(2),imgsizePerm(3),imgsizePerm(4),imgsizePerm(5));
                    elseif(length(imgsize) == 6)
                        img = reshape(permute(img,[2 1 3 4]),imgsizePerm(1),imgsizePerm(2),imgsizePerm(3),imgsizePerm(4),imgsizePerm(5),imgsizePerm(6));
                    end
                    img = ipermute(img, para.permVecTrafo(i,:));
                    
%                     % t-y-z-x-cha
%                     img = permute(img,[2 3 1 4 5]); % y-z-t-x-cha
%                     tmp = permute(reshape(img,[imgsize(2)*imgsize(3), imgsize(1), imgsize(4), imgsize(5)]),[2 1 3 4]);
%                     tmp = permute(mtimesx(para.pc{i},tmp),[2 1 3 4]); % t - y*z - x - cha
%                     img = permute(reshape(tmp,[imgsize(2), imgsize(3), imgsize(1), imgsize(4), imgsize(5)]),[3 1 2 4 5]);
                    
%                     if(all(size(img) == para.size))
%                         img = ipermute(mtimesx(para.pc{i},permute(img,para.permVecTrafo(i,1:n_dim))),para.permVecTrafo(i,1:n_dim));
%                     else
%                         img = ipermute(mtimesx(para.pcKernel{i},permute(img,para.permVecTrafo(i,1:n_dim))),para.permVecTrafo(i,1:n_dim));
%                     end
                end
                
            case 'wavelet_lab'
                % wavelet => cart
                imgsize = ones(1,5);
                imgTmp = size(img);
                imgsize(1:length(imgTmp)) = imgTmp;
                for i=1:imgsize(5)
                    for j=1:imgsize(4)
                        for k=1:imgsize(3)
                            img(:,:,k,j,i) = para.W' * img(:,:,k,j,i);
                        end
                    end
                end
                if(para.zeroPad)
                    img = crop(img,[para.size]);
                else
                    imgOut = zeros(para.size,para.precision);
                    for i=1:imgsize(5)
                        for j=1:imgsize(4)
                            for k=1:imgsize(3)
                                imgOut(:,:,k,j,i) = imresize(img(:,:,k,j,i), [para.size(1) para.size(2)], para.rescaleInterp);
                            end
                        end
                    end
                    img = imgOut;
                    clear 'imgOut';
                end
                
            case 'wavelet_mat'
                % wavelet => cart
                imgIn = getMeta(img);
                data = getData(img);
                img = zeros(para.size,para.precision);
                imgsize = para.imgsize;
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:para.size(5) % cha
                        for i=1:para.size(4)
                            if(para.shape == 3)
                                imgIn{h}.dec = data{h};
                                img(:,:,:,i,h) = waverec3(imgIn{h});
                                dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                                continue;
                            end
                            for j=1:para.size(3) % freq
                                if(para.shape == 2)
                                    img(:,:,j,i,h) = waverec2(data{h}(:,:,j,i), imgIn{h}(:,:,j,i), para.wFilter);
                                    dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                                    continue;
                                end
                                for k=1:para.size(2) % z
                                    if(para.shape == 1)
                                        img(:,k,j,i,h) = waverec(data{h}(:,k,j,i), imgIn{h}(:,k,j,i), para.wFilter);
                                    end
                                    dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + (j-1)*size(img,2) + k)/prod(imgsize(2:end)));
                                end
                                dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                            end
                            dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                        end
                        dispProgress('Trafo', (h*size(img,4)*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:para.size(4) % cha                      
                        if(para.shape == 3)   
                            imgIn{i}.dec = data{i};
                            img(:,:,:,i) = waverec3(imgIn{i});
                            dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                            continue;
                        end
                        for j=1:para.size(3) % freq
                            if(para.shape == 2)
                                img(:,:,j,i) = waverec2(data{i}(:,:,j), imgIn{i}(:,:,j), para.wFilter);
                                dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                                continue;
                            end
                            for k=1:para.size(2) % z
                                if(para.shape == 1)
                                    img(:,k,j,i) = waverec(data{i}(:,k,j), imgIn{i}(:,k,j), para.wFilter);
                                end
                                dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + (j-1)*size(img,2) + k)/prod(imgsize(2:end)));
                            end
                            dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                        end
                        dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:para.size(3)
                        if(para.shape == 2)
                            img(:,:,i) = waverec2(data{i}(:,:), imgIn{i}(:,:), para.wFilter);
                            dispProgress('Trafo', (i*size(img,2))/prod(imgsize(2:end)));
                            continue;
                        end
                        for j=1:para.size(2)
                            if(para.shape == 1)
                                img(:,j,i) = waverec(data{i}(:,j), imgIn{i}(:,j), para.wFilter);
                                dispProgress('Trafo', ((i-1)*size(img,2) + j)/prod(imgsize(2:end)));
                            end
                        end
                        dispProgress('Trafo', (i*size(img,2))/prod(imgsize(2:end)));
                    end
                end
                dispProgress('Trafo', 'Close');
                
            case 'mellin'
                % mellin => cart      
                if(para.shape == 2)                    
                    imgIn = getData(img); %2D trafo allows arbitrary transformation size
                    imgsize = getMeta(img);
                else
                    imgIn = img;
                    imgsize = para.imgsize;
                end
                img = zeros(para.size,para.precision);
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:size(img,5) % cha
                        for i=1:size(img,4)
                            if(para.shape == 3)
                                img(:,:,:,i,h) = AFMT3(imgIn(:,:,:,i,h), para.sigma, size(imgIn,1), size(imgIn,2), size(imgIn,3), para.extrapVal, para.interp, 'iF-AFMT');
                                dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                                continue;
                            end
                            for j=1:size(img,3) % freq
                                if(para.shape == 2)
                                    img(:,:,j,i,h) = AFMT2(imgIn{h}(:,:,j,i), para.sigma, imgsize(1), imgsize(2), para.extrapVal, para.interp, 'iF-AFMT', para.bTrafoType);
                                    dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                                end
                            end
                            dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                        end
                        dispProgress('Trafo', (h*size(img,4)*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:size(img,4) % cha
                        if(para.shape == 3)
                            img(:,:,:,i) = AFMT3(imgIn(:,:,:,i), para.sigma, size(imgIn,1), size(imgIn,2), size(imgIn,3), para.extrapVal, para.interp, 'iF-AFMT');
                            dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                            continue;
                        end
                        for j=1:size(img,3) % freq
                            if(para.shape == 2)   
                                img(:,:,j,i) = AFMT2(imgIn{i}(:,:,j), para.sigma, imgsize(1), imgsize(2), para.extrapVal, para.interp, 'iF-AFMT', para.bTrafoType);
                                dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                            end
                        end
                        dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:size(img,3)
                        if(para.shape == 2)
                            img(:,:,i) = AFMT2(imgIn{i}(:,:), para.sigma, imgsize(1), imgsize(2), para.extrapVal, para.interp, 'iF-AFMT', para.bTrafoType);
                        end
                        dispProgress('Trafo', (i*size(img,2))/prod(imgsize(2:end)));
                    end
                end
                dispProgress('Trafo', 'Close');
                
            case 'curvelet'
                imgIn = getData(img);
                imgsize = getMeta(img);                
                img = zeros(imgSize,para.precision); %initialize backward transformed image
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:imgsize(5) % cha
                        for i=1:imgsize(4)
                            if(para.shape == 3)
                                img(:,:,:,i,h) =  fdct3d_inverse_mex(imgsize(1), imgsize(2), imgsize(3), para.nbscales, para.nbdstz_coarse, para.allCurvelets, imgIn{h}{i});                
                                dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                                continue;
                            end
                            for j=1:imgsize(3) % freq
                                if(para.shape == 2)
                                    img(:,:,j,i,h) = ifdct_usfft_mex(imgsize(1), imgsize(2), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,j,i,h));
                                    dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                                end
                            end
                            dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                        end
                        dispProgress('Trafo', (h*size(img,4)*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:imgsize(4) % cha
                        if(para.shape == 3)
                            img(:,:,:,i) =  fdct3d_inverse_mex(imgsize(1), imgsize(2), imgsize(3), para.nbscales, para.nbdstz_coarse, para.allCurvelets, imgIn{i});                
                            dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                            continue;
                        end
                        for j=1:imgSize(3)
                            if(para.shape == 2)
                                img(:,:,j,i) = ifdct_usfft_mex(imgsize(1), imgsize(2), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,j,i));
                                dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                            end
                        end
                        dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:size(img,3)
                        if(para.shape == 2)
                            img(:,:,i) = ifdct_usfft_mex(imgsize(1), imgsize(2), para.nbscales, para.nbdstz_coarse, para.allCurvelets, img(:,:,i));
                            dispProgress('Trafo', (i*size(img,2))/prod(imgsize(2:end)));
                        end
                    end
                end
                dispProgress('Trafo', 'Close');
                                            
            case 'surfacelet'
                % surfacelet => cart
                imgIn = getData(img);
                meta = getMeta(img);
                recinfo = meta.recinfo;
                idx = recinfo{1,2};
                imgsize = para.imgsize;
                clear 'img'
                img = zeros([para.padsize*ones(1,para.shape),para.imgsize(para.shape+1:end)],para.precision);
                dispProgress('Trafo', 0, prod(imgsize(2:end)));
                if(strcmp(para.dimensionality,'4D')) % t-y-z-x-cha
                    for h=1:size(img,5) % cha
                        for i=1:size(img,4)
                            if(para.shape == 3)
%                                 [tmp, idx] = flattenCellMatrix(imgIn{1,h}{1,i});
                                imgReal = reconFlatCellMatrix(cellfun(@(x) real(x), imgIn{1,h}{1,i}, 'UniformOutput', false),idx);
                                imgImag = reconFlatCellMatrix(cellfun(@(x) imag(x), imgIn{1,h}{1,i}, 'UniformOutput', false),idx);
                                clear 'tmp';
                                img(:,:,:,i,h) = Surfrec(imgReal,recinfo{1,1}) + 1i * Surfrec(imgImag,recinfo{1,1});
                                clear 'imgReal' 'imgImag'
                                dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                                continue;                                
                            end
                            for j=1:size(img,3)
                                if(para.shape == 2)
%                                     [tmp, idx] = flattenCellMatrix(imgIn{1,h}{1,i}(:,j));
                                    imgReal = reconFlatCellMatrix(cellfun(@(x) real(x), imgIn{1,h}{1,i}(:,:,j), 'UniformOutput', false),idx);
                                    imgImag = reconFlatCellMatrix(cellfun(@(x) imag(x), imgIn{1,h}{1,i}(:,:,j), 'UniformOutput', false),idx);
                                    clear 'tmp' 'idx'
                                    img(:,:,j,i,h) = Surfrec(imgReal,recinfo{1,1}) + 1i * Surfrec(imgImag,recinfo{1,1});
                                    clear 'imgReal' 'imgImag'
                                    dispProgress('Trafo', ((h-1)*size(img,4)*size(img,3)*size(img,2) + (i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                                end
                            end
                        end
                    end
                elseif(strcmp(para.dimensionality,'3D') || strcmp(para.dimensionality,'2Dt'))
                    for i=1:size(img,4) % cha
                        if(para.shape == 3)
%                             [tmp, idx] = flattenCellMatrix(imgIn{1,i});
                            imgReal = reconFlatCellMatrix(cellfun(@(x) real(x), imgIn{1,i}, 'UniformOutput', false),idx);
                            imgImag = reconFlatCellMatrix(cellfun(@(x) imag(x), imgIn{1,i}, 'UniformOutput', false),idx);
                            clear 'tmp';
                            img(:,:,:,i) = Surfrec(imgReal,recinfo{1,1}) + 1i * Surfrec(imgImag,recinfo{1,1});
                            clear 'imgReal' 'imgImag'
                            dispProgress('Trafo', (i*size(img,3)*size(img,2))/prod(imgsize(2:end)));
                            continue;
                        end
                        for j=1:size(img,3)
                            if(para.shape == 2)
%                                 [tmp, idx] = flattenCellMatrix(imgIn{1,i}(:,j));
                                imgReal = reconFlatCellMatrix(cellfun(@(x) real(x), imgIn{1,i}(:,:,j), 'UniformOutput', false),idx);
                                imgImag = reconFlatCellMatrix(cellfun(@(x) imag(x), imgIn{1,i}(:,:,j), 'UniformOutput', false),idx);
                                clear 'tmp';
                                img(:,:,j,i) = Surfrec(imgReal,recinfo{1,1}) + 1i * Surfrec(imgImag,recinfo{1,1});
                                clear 'imgReal' 'imgImag'
                                dispProgress('Trafo', ((i-1)*size(img,3)*size(img,2) + j*size(img,2))/prod(imgsize(2:end)));
                            end
                        end
                    end
                elseif(strcmp(para.dimensionality,'2D'))
                    for i=1:size(img,3) % cha
                        if(para.shape == 2)
%                             [tmp, idx] = flattenCellMatrix(imgIn{1,i});
                            imgReal = reconFlatCellMatrix(cellfun(@(x) real(x), imgIn{1,i}, 'UniformOutput', false),idx);
                            imgImag = reconFlatCellMatrix(cellfun(@(x) imag(x), imgIn{1,i}, 'UniformOutput', false),idx);
                            clear 'tmp';
                            img(:,:,i) = Surfrec(imgReal,recinfo{1,1}) + 1i * Surfrec(imgImag,recinfo{1,1});
                            clear 'imgReal' 'imgImag'
                            dispProgress('Trafo', (i*size(img,2))/prod(imgsize(2:end)));
                        end
                    end
                end
                dispProgress('Trafo', 'Close');
                clear 'imgIn'
                if(para.zeroPad)
                    img = crop(img,para.imgsize);
                else
                    imgOut = zeros(para.imgsize,para.precision);
                    if(para.shape == 3)
                        for i=1:size(img,5)
                            for j=1:size(img,4)
                                RI = imref3d(size(img(:,:,:,j,i)),1,1,1);
                                scaleFactor = para.size(1:3)./size(img(:,:,:,j,i));
                                tFormResize = affine3d([scaleFactor(2) 0 0 0; 0 scaleFactor(1) 0 0; 0 0 scaleFactor(3) 0; 0 0 0 1]);
                                imgOut(:,:,:,j,i) = crop(imwarp(img(:,:,:,j,i),RI,tFormResize,para.rescaleInterp),para.size(1:3));
                            end
                        end
                    elseif(para.shape == 2)
                        for i=1:size(img,5)
                            for j=1:size(img,4)
                                for k=1:size(img,3)
                                    imgOut(:,:,k,j,i) = imresize(img(:,:,k,j,i), para.imgsize(1:2), para.rescaleInterp);
                                end
                            end
                        end
                    end
                    img = imgOut;
                    clear 'imgOut' 'idx';
                end
        end
        
        if(~isempty(para.permRule))
            img = ipermute(img,para.permRule);
        end
        if(para.kspaceTrafo && ~isempty(para.fftdim))
            img = fftnshift(img,para.fftdim,para.scrambledim);
        end
    end
end
