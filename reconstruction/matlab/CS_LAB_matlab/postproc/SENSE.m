function imgSENSE = SENSE(obj, image, postproc)
%SENSE Summary of this function goes here
%  image:   cell array 1 x nCha
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

%% extract dimensions and prepare outputs
[nPha, nFreq, nZ, nTime, nCha] = obj.measPara.dim;

if(strcmp(obj.measPara.dimension,'2D'))
    if(nTime == 1)
        % y-x-cha
        imgSENSE = complex(zeros(nPha,nFreq),zeros(nPha,nFreq));
        imgIn = cell2mat(shiftdim(image,-1)); % y-x-cha     
        kSpace = imgIn;
    else
        % y-x-t-cha    
        imgSENSE = complex(zeros(nPha,nFreq,nTime),zeros(nPha,nFreq,nTime));
        imgIn = cell2mat(shiftdim(image,-2)); % y-x-t-cha
        kSpace = permute(imgIn,[1 2 4 3]); % y-x-cha-t
    end
    
elseif(strcmp(obj.measPara.dimension,'3D'))
    imgSENSE = complex(zeros(nPha,nFreq,nZ),zeros(nPha,nFreq,nZ));
    imgIn = cell2mat(shiftdim(image,-2)); % y-x-z-cha
    
elseif(strcmp(obj.measPara.dimension,'4D'))
    imgSENSE = complex(zeros(nPha,nFreq,nZ,nTime),zeros(nPha,nFreq,nZ,nTime));
    imgIn = cell2mat(shiftdim(image,-3)); % y-x-z-t-cha
    
else
    error('SENSE(): Undetermined dimension');
end      
  

%% create calibration matrix from full k-space and
%% calculate SENSE maps and correct images
if(strcmp(obj.measPara.dimension,'2D'))      
    kSpace = fft2shift(kSpace);
    kernelSize = obj.kernelSize;

    dispProgress('SENSE - Time', 0, nTime);
    for t=1:nTime
        if(isempty('obj.calibSize','var'))
            % shouldn't happen, but to be absolutely sure
            calibSize = obj.calibrationSize(obj.fullMask(:,:,t)); % k_y - k_x
        else
            if(iscell(obj.calibSize))
                calibSize = obj.calibSize{obj.currSlice,t};
            else
                calibSize = obj.calibSize;
            end
        end

        kCalib = crop(kSpace(:,:,:,t), [calibSize, nCha]);

        A = calibMatrix2D(kCalib,kernelSize,nCha); 
        
        senseMap = calcSENSEMap2D(A,nPha,nFreq,nCha,kernelSize,postproc);

        if(nTime == 1)
            imgSENSE = sum(senseMap.*imgIn,3);
        else
            imgSENSE(:,:,t) = sum(senseMap.*permute(imgIn(:,:,t,:),[1 2 4 3]),3);
        end 
        dispProgress('SENSE - Time', t/nTime);
    end
    dispProgress('SENSE - Time', 'Close');

elseif(strcmp(obj.measPara.dimension,'3D'))        
% %         kSpace = permute(imgIn, [1 3 2 4]); % y-z-x-cha
    kSpace = fftnshift(imgIn); % k_y-k_x-k_z-cha

% %         kernelSize = [obj.kernelSize(1), obj.kernelSize(3)];
    calibSize = obj.calibSize;

% %         kCalib = crop(kSpace, [calibSize, nCha]);

% %         A = cell(nZ, 1);
% %         for z=1:nZ
% %             kCalib = crop(permute(kSpace(:,:,z,:),[1 2 4 3]), [calibSize(1), calibSize(3), nCha]);
% %             A{z,1} = calibMatrix2D(kCalib,kernelSize,nCha);
% %         end
    kernelSize = [obj.kernelSize(1) obj.kernelSize(3) obj.kernelSize(2)];
    kCalib = crop(kSpace, [calibSize(1), calibSize(3), calibSize(2), nCha]);
    A = calibMatrix3D(kCalib,kernelSize,nCha);
    
    senseMap = calcSENSEMap3D(A,nPha,nFreq,nZ,nCha,kernelSize,postproc);
        
    imgSENSE = sum(senseMap.*imgIn,4);

elseif(strcmp(obj.measPara.dimension,'4D'))        
    kSpace = permute(imgIn, [1 2 3 5 4]); % y-x-z-cha-t
    kSpace = fftnshift(kSpace,1:3);
    kernelSize = [obj.kernelSize(1) obj.kernelSize(3) obj.kernelSize(2)];

    dispProgress('SENSE - Time', 0, nTime);
    for t=1:nTime
        if(iscell(obj.calibSize))
            calibSize = obj.calibSize{t};
        else
            calibSize = obj.calibSize;
        end

        kCalib = crop(kSpace(:,:,:,:,t), [calibSize(1), calibSize(3), calibSize(2), nCha]);
        A = calibMatrix3D(kCalib,kernelSize,nCha);
        
        senseMap = calcSENSEMap3D(A,nPha,nFreq,nZ,nCha,kernelSize,postproc);
        
        imgSENSE(:,:,:,t) = sum(senseMap.*permute(imgIn(:,:,:,t,:),[1 2 3 5 4]),4);
% % %             for z=1:nTime %%??????????????????
% % %                 A{z,t} = calibMatrix3D(kCalib,kernelSize,nCha); 
% % %             end
        dispProgress('SENSE - Time', t/nTime);
    end
    dispProgress('SENSE - Time', 'Close');
    
else
    error('SENSE(): Undetermined dimension');
end 
    
end





function A = calibMatrix2D(kCalib,kernelSize,nCha)
% calibration Matrix 2D
A = [];

for n=1:nCha
	tmp  = im2col(kCalib(:,:,n),kernelSize,'sliding').';
	A = [A, tmp];
end
end

function A = calibMatrix3D(kCalib,kernelSize,nCha)
% calibration Matrix 3D
A = [];

for n=1:nCha
    if(isreal(kCalib(:,:,:,n)))
        tmp = im3colR(kCalib(:,:,:,n),kernelSize).';
    else
        tmp = im3colC(kCalib(:,:,:,n),kernelSize).';
    end
	A = [A, tmp];
end
end


function senseMap = calcSENSEMap2D(A,nPha,nFreq,nCha,kernelSize,postproc)
% calculate the senseMap for 2D case

global prop

openPool = prop.openPool;

[~,S,V] = svd(A);

% find parallel part of V -> determine rank of A
%         tol = max(size(A))*eps(max(S(:)));
rA = sum(sum(S > postproc.tol)); % ^= K
% rA = 80;
Vpara = V(:,1:rA);
% kfilter = cell(1,rA);
cols = cell(1,rA);
Gr = cell(nPha,nFreq);

setMatlabPool('open');
dispProgress('SENSE - Kernel', 0, rA);

if(openPool)
    parfor j=1:rA
        cols{j} = calcSENSEMap2DKernel(j,Vpara,kernelSize,nPha,nFreq,nCha);
        dispProgress('SENSE - Kernel', j);
    end
    setMatlabPool('close');
else
    for j=1:rA
        cols{j} = calcSENSEMap2DKernel(j,Vpara,kernelSize,nPha,nFreq,nCha);
        dispProgress('SENSE - Kernel', j/rA);
    end
end
dispProgress('SENSE - Kernel', 'Close');

dispProgress('SENSE - Concatenate', 0, rA);
for j=1:rA
    Gr = cellfun(@(x,y) [x;y], Gr, cols{j}, 'UniformOutput', false); 
    dispProgress('SENSE - Concatenate', j/rA);
end
dispProgress('SENSE - Concatenate', 'Close');

% final G = Gr^H * Gr
G = cellfun(@(x) x'*x, Gr, 'UniformOutput', false);
[Veig, ~] = cellfun(@(x) eig(x,'nobalance'), G, 'UniformOutput', false);

eigMaps = cellfun(@(x) x(:,nCha), Veig, 'UniformOutput', false);
eigMaps = cell2mat(cellfun(@(x) shiftdim(x,-2), eigMaps, 'UniformOutput', false));

% create reference
phases = angle(eigMaps);

eigMaps = abs(eigMaps).*exp(1i*(phases-repmat(phases(:,:,nCha),[1,1,nCha])));

fctr = sqrt(sum(conj(eigMaps).*eigMaps,3)).^(-1);
senseMapHelper = repmat(fctr,[1 1 nCha]).*eigMaps;
%         senseMap{z,t} = senseMapHelper.*conj(senseMapHelper);
senseMap = conj(senseMapHelper);

end


function cols = calcSENSEMap2DKernel(j,Vpara,kernelSize,nPha,nFreq,nCha)
%     fprintf('Loop %d/%d\n',j,rA);
%     helper = permute(reshape(Vpara(:,j),[kernelSize,nCha]), [2 1 3]);
helper = reshape(Vpara(:,j),[kernelSize,nCha]);
%             kfilter{j} = crop(ifft2shift(zpad(helper, nPha+obj.kernelSize(1)-1, nFreq+obj.kernelSize(2)-1, nCha)),nPha,nFreq,nCha);
kfilter = ifft2shift(zpad(helper(end:-1:1,end:-1:1,:), nPha, nFreq, nCha));
%             kfilter{j} = GFFT' * zpad(helper, nPha, nFreq, nCha);
helper = mat2cell(kfilter, ones(1,nPha), ones(1,nFreq), nCha);
cols = cellfun(@(x) permute(x,[1 3 2]), helper, 'UniformOutput', false);

end


function senseMap = calcSENSEMap3D(A,nPha,nFreq,nZ,nCha,kernelSize,postproc)
% calculate the senseMap for 3D case

global prop

openPool = prop.openPool;

[~,S,V] = pcaSVD(A,size(A,2));

% find parallel part of V
% -> determine rank of A
%         tol = max(size(A))*eps(max(S(:)));
helper = S./max(S(:));
rA = sum(sum(helper > postproc.tol)); % ^= K
% rA = 80;
Vpara = V(:,1:rA);
%         kfilter = cell(1,rA);
cols = cell(1,rA);
Gr = cell(nPha,nFreq,nZ);

setMatlabPool('open');
dispProgress('SENSE - Kernel', 0, rA);

if(openPool)
    parfor j=1:rA
        cols{j} = calcSENSEMap3DKernel(j,Vpara,kernelSize,nPha,nFreq,nZ,nCha);
        dispProgress('SENSE - Kernel', j);
    end
    setMatlabPool('close');
else
    for j=1:rA
        cols{j} = calcSENSEMap3DKernel(j,Vpara,kernelSize,nPha,nFreq,nZ,nCha);
        dispProgress('SENSE - Kernel', j/rA);
    end
end
dispProgress('SENSE - Kernel', 'Close');

dispProgress('SENSE - Concatenate', 0, rA);
for j=1:rA
    Gr = cellfun(@(x,y) [x;y], Gr, cols{j}, 'UniformOutput', false); 
    dispProgress('SENSE - Concatenate', j/rA);
end
dispProgress('SENSE - Concatenate', 'Close');

% final G = Gr^H * Gr
G = cellfun(@(x) x'*x, Gr, 'UniformOutput', false);
[Veig, ~] = cellfun(@(x) eig(x,'nobalance'), G, 'UniformOutput', false);

eigMaps = cellfun(@(x) x(:,nCha), Veig, 'UniformOutput', false);
eigMaps = cell2mat(cellfun(@(x) shiftdim(x,-3), eigMaps, 'UniformOutput', false));

% create reference
phases = angle(eigMaps);

eigMaps = abs(eigMaps).*exp(1i*(phases-repmat(phases(:,:,:,nCha),[1,1,1,nCha])));

fctr = sqrt(sum(conj(eigMaps).*eigMaps,4)).^(-1);
senseMapHelper = repmat(fctr,[1 1 1 nCha]).*eigMaps;
%         senseMap{z,t} = senseMapHelper.*conj(senseMapHelper);
senseMap = conj(senseMapHelper);

end

function cols = calcSENSEMap3DKernel(j,Vpara,kernelSize,nPha,nFreq,nZ,nCha)
helper = reshape(Vpara(:,j),[kernelSize,nCha]);
%             kfilter{j} = crop(ifft2shift(zpad(helper, nPha+obj.kernelSize(1)-1, nFreq+obj.kernelSize(2)-1, nCha)),nPha,nFreq,nCha);
kfilter = ifftnshift(zpad(helper(end:-1:1,end:-1:1,end:-1:1,:), nPha, nFreq, nZ, nCha),1:3);
%             kfilter{j} = GFFT' * zpad(helper, nPha, nFreq, nCha);
helper = mat2cell(kfilter, ones(1,nPha), ones(1,nFreq), ones(1,nZ), nCha);
cols = cellfun(@(x) permute(x,[1 4 2 3]), helper, 'UniformOutput', false);
end