function [imageCha] = recon_main(obj, iRep, iAvg)
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
% nSlices = obj.measPara.dim(3); % size(obj.kSpace,1);
% nSlices = obj.measPara.LCall(1);
% nTime = obj.measPara.dim(4);
nSlices = size(obj.kSpace,1);
imageCha = cell(nSlices,nCha);

% prepare wavelet filter names
waveletFilter_l1 = obj.trafo.waveletFilter;
waveletFilterSize_l1 = obj.trafo.waveletFilterSize;
waveletFilter_l12 = obj.trafo.waveletFilter_l12;
waveletFilterSize_l12 = obj.trafo.waveletFilterSize;
if(strcmp(waveletFilter_l1,'dmey'))
    obj.trafo.waveletFilterName_l1 = waveletFilter_l1;            
else
    obj.trafo.waveletFilterName_l1 = [waveletFilter_l1,num2str(waveletFilterSize_l1)];
end
if(strcmp(waveletFilter_l12,'dmey'))
    obj.trafo.waveletFilterName_l12 = waveletFilter_l12;            
else
    obj.trafo.waveletFilterName_l12 = [waveletFilter_l12,num2str(waveletFilterSize_l12)];
end
      
if(strcmp(obj.measPara.dimension,'2D')) % atm not really needed, since all test spaces are 3D.
    % 2D case: x-y-t or x-y (sparse dimensions: y-t)
    	
    dispProgress('Slices', 0, nSlices);
    for iSli = 1:nSlices      
        % prepare k-space 
        % => k_y - k_x - t - sli - cha
        kSpaceL = squeeze(cell2mat(shiftdim(obj.kSpace(iSli,:,iRep,iAvg),-3)));
%         if((obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1 && nSlices == 1) || ...
%                 (iRep == size(obj.kSpace,3) && iAvg == size(obj.kSpace,4) && iSli == nSlices) || ...
%                 (obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1 && iSli == nSlices))
%             obj.kSpace = [];
%         end
        if(strcmp(obj.measPara.precision,'double'))
            kSpaceL = double(kSpaceL);
        end
        % extraction of sampling mask
        % !!! Must be the first operation, because fft and ifft operations,
        % destroy this information !!!
        % same undersampling for all coils and fully sampled in k_x direction
        % => reduce 4D problem (k_y - k_x - t - cha) to 2D problem (k_y - t)
    %     lMask = cellfun(@(x) abs(x) > 0, all_KSpaces, 'UniformOutput', false);
    %     obj.fullMask = lMask{1}; % k_y - k_x - t 
        obj.fullMask = abs(kSpaceL) > 0; % k_y - k_x - t - sli - cha
        if(obj.measPara.oversampling{2,1})
            obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :, :); % anti-aliasing
        %     mask = abs(sum(permute(fullMask,[3 1 2]),3)) > 0;

            kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - t - cha
            kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :, :); % anti-aliasing
            kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - t -cha
            kSpaceL = kSpaceL .* obj.fullMask; % just take acquired points
            obj.measPara.dim(2) = size(kSpaceL,2);
        end
        
        % [nPha, nFreq, nTime, nCha] = size(kSpaceL);
        
        %  without change: k_y - k_x - k_z - cha  
%         % prepare k-space and mask => k_y - k_x - cha - t
        kSpaceL = permute(kSpaceL, [1 2 4 3]);
        obj.fullMask = permute(obj.fullMask, [1 2 4 3]);

        % CS 2D reconstruction
        fprintf('CS 2D reconstruction\n');
        imgOut = obj.CS_2D(kSpaceL);
        
        imgOut = cell2mat(shiftdim(imgOut,-2)); % y-x-sli-cha (t not supported yet)
        imageCha(iSli,:) = shiftdim(mat2cell(imgOut,size(imgOut,1),size(imgOut,2),ones(1,size(imgOut,3)), ones(1,size(imgOut,4))),2);
        dispProgress('Slices', iSli/nSlices);
    end
    dispProgress('Slices','Close');
    imageCha = cellfun(@(x) ipermute(x,[1 2 4 3]), imageCha, 'UniformOutput', false);
    imageCha = cellfun(@(x) flipdim(x,1), imageCha, 'UniformOutput', false);
    
elseif(strcmp(obj.measPara.dimension,'3D'))
    % 3D case: x-y-z (sparse dimensions: y-z)
           
    % prepare k-space 
    % => k_y - k_x - k_z - cha
    kSpaceL = cell2mat(shiftdim(obj.kSpace(:,:,iRep,iAvg),-2));
    if((obj.measPara.LCall(4) == 1 && obj.measPara.LCall(2) == 1) || (iRep == size(obj.kSpace,3) && iAvg == size(obj.kSpace,4)))
        obj.kSpace = [];
    end
    if(strcmp(obj.measPara.precision,'double'))
        kSpaceL = double(kSpaceL);
    end
    
    % extraction of sampling mask
    % !!! Must be the first operation, because fft and ifft operations,
    % destroy this information !!!
    % same subsampling for all channels
    obj.fullMask = abs(kSpaceL) > 0; % -> k_y-k_x-k_z-cha
    
    if(length(unique(sum(obj.fullMask,2))) > 2) % ensure same subsampling along readout direction
        val = max(squeeze(sum(sum(sum(obj.fullMask,1),2),3)))/size(obj.fullMask,2);
        i = 1;
        while(i <= size(obj.fullMask,2) && sum(sum(obj.fullMask(:,i,:,1))) ~= val)
            i = i+1;
        end
        if(i > size(obj.fullMask,2))
            error('recon_main(): Incomplete undersampling mask!');
        end
        obj.fullMask = repmat(obj.fullMask(:,i,:,:), [1 size(obj.fullMask,2) 1 1]);
    end
        
    if(obj.measPara.oversampling{2,1})
        obj.fullMask = obj.fullMask(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
        kSpaceL = fftshift(ifft(ifftshift(kSpaceL, 2),[],2),2); % -> k_y - x - k_z - cha
        kSpaceL = kSpaceL(:, obj.measPara.oversampling{1,1}, :, :); % anti-aliasing
        kSpaceL = fftshift(fft(ifftshift(kSpaceL,2),[],2),2); % -> k_y - k_x - k_z - cha
        kSpaceL = kSpaceL .* obj.fullMask; % just take acquired points
        obj.measPara.dim(2) = size(kSpaceL,2);
    end
    
    % remove all unnecessary variables
    clear 'kSpace' 'kernel' 'window' 'AtA';
    kSpaceL = permute(kSpaceL,[1 3 2 4]); 
    obj.fullMask = permute(obj.fullMask,[1 3 2 4]);
    
    % CS 3D reconstruction
    if strcmp(obj.reconDIM,'2D')
        fprintf('CS 2D reconstruction of 3D data\n');
        imageCha = obj.CS_2D(kSpaceL);
    elseif strcmp(obj.reconDIM,'3D')
        fprintf('CS 3D reconstruction of 3D data\n');
        imageCha = obj.CS_3D(kSpaceL);
    end;
    imageCha = cellfun(@(x) ipermute(x,[1 3 2 4]), imageCha, 'UniformOutput', false);
    imageCha = cellfun(@(x) flipdim(flipdim(x,1),2), imageCha, 'UniformOutput', false);
    
elseif(strcmp(obj.measPara.dimension,'4D'))
    % to do
elseif(strcmp(obj.measPara.dimension,'5D'))
    % to do
end;

end