function imageCha = recon_main(obj, iRep, iAvg)
% main reconstruction function
% prepare constraints and select algorithm depending on image dimensionality
%
% input:
% obj       CS reconstruction object (holding all parameters)
% iRep      current repetition
% iAvg      current average
%
% output:
% imageCha  individual reconstructed channel images
%
% (c) Thomas Kuestner 
% -------------------------------------------------------------------------

nCha = obj.measPara.dim(5);
nSlices = size(obj.kSpace,1);
nTime = obj.measPara.dim(4);
imageCha = cell(nSlices,nCha);

if(strcmp(obj.measPara.dimension,'2D'))

    dispProgress('Slices', 0, nSlices);
    for iSli = 1:nSlices
        % prepare k-space 
        % => k_y - k_x - t - cha
        kSpaceL = cell2mat(shiftdim(obj.kSpace(iSli,:,iRep,iAvg),-2));
        if(strcmp(obj.measPara.precision,'double'))
            kSpaceL = double(kSpaceL);
        end
        % same undersampling for all coils and fully sampled in k_x direction
        % => reduce 4D problem (k_y - k_x - t - cha) to 2D problem (k_y - t)
        obj.fullMask = abs(kSpaceL(:,:,:,1)) > 0; % k_y - k_x - t
        if(obj.measPara.oversampling{2,1})
            kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - t - cha
            kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
            obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :);
            kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - t -cha
            kSpaceL = kSpaceL .* repmat(obj.fullMask,[1 1 1 size(kSpaceL,4)]);
            obj.measPara.dim(2) = size(kSpaceL,2);
        end
        [nPha, nFreq, nTime, nCha] = size(kSpaceL); 
   
        % ensure that empty entries are zero
        mask = false(nPha,nFreq,nCha,nTime);
        for i=1:nTime
            mask(:,:,:,i) = repmat(obj.fullMask(:,:,i),[1,1,nCha]); 
        end
        
        imageCha(iSli,:) = shiftdim(mat2cell(ifftnshift(kSpaceL,1:2,1:2),nPha,nFreq,nTime,ones(1,nCha)),2);
        
        dispProgress('Slices', iSli/nSlices);
    end
    dispProgress('Slices', 'Close');

elseif(strcmp(obj.measPara.dimension,'3D'))
        
    % prepare k-space 
    % => k_y - k_z - k_x - cha
    kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-2));
    kSpaceL = permute(kSpaceL,[1 3 2 4]);
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    % same undersampling for all coils and fully sampled in k_x direction
    obj.fullMask = abs(kSpaceL) > 0; % k_y - k_z - k_x - cha
    if(obj.measPara.oversampling{2,1})
        kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 3),[],3),3); % -> k_y - k_z - x - cha
        kSpaceL = kSpaceL(:, :, obj.measPara.oversampling{1,1}, :); % anti-aliasing
        obj.fullMask = obj.fullMask(:,:, obj.measPara.oversampling{1,1}, :);
        kSpaceL = fftshift(fft(ifftshift(kSpaceL,3),[],3),3); % -> k_y - k_z - k_x - cha
        kSpaceL = kSpaceL .* obj.fullMask;
        obj.measPara.dim(2) = size(kSpaceL,3);
    end
    
    [nPha, nZ, nFreq, nCha] = size(kSpaceL);   
    kSpaceL = permute(kSpaceL,[1 3 2 4]); % y-x-z-cha
           
	imageCha = shiftdim(mat2cell(ifftnshift(kSpaceL,1:3,1:3), nPha, nFreq, nZ, ones(1,nCha)),2);        
            
elseif(strcmp(obj.measPara.dimension,'4D'))
    % TODO

end
end