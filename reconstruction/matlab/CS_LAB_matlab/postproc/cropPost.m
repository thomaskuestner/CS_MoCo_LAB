function image = cropPost(image,oversampling,aniso,interp,turnImage,FreqOversamplingCorr,measPara)
%CROPPOST crop and flip image
% 2D: y-x-slices
% 2Dt:y-x-t-slices
% 3D: y-x-z
% 4D: y-x-z-t
% 5D: y-x-z-t-g
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% 1.) correct anisotropy
% kspacesamples -> ftlen (just needed for phase encoding directions (z and
% y)
if(isempty(aniso))
    aniso = size(image);
end
% if(aniso(1) > aniso(2)) % anisotropic y resolution never gets larger than readout (x)
%     aniso(1) = aniso(2);
% end
% tmp = aniso(~ismember(aniso,1));
if(any(aniso(1:ndims(image)) ~= size(image)))
    imageOut = zeros(aniso(1),aniso(3),aniso(2),aniso(4),aniso(5));
    image = permute(image,[1 3 2 4 5]);
    for i=1:size(image,3)
        for j=1:size(image,4)
            imageOut(:,:,i,j,:) = imresize(image(:,:,i,j,:), [aniso(1), aniso(3)], interp);
        end
    end
    image = ipermute(imageOut, [1 3 2 4 5]);
end
    
% correct oversampling
if(FreqOversamplingCorr)
    oversampling{2,1} = ~oversampling{2,1}; % always correct readout oversampling
end
oversamplingLoc = cell(1,3);
tmpPerm = [2 1 3];
for i=1:3
    if(~oversampling{2,i})
        oversamplingLoc{1,i} = 1:size(image,tmpPerm(i));
    else
        oversamplingLoc{1,i} = oversampling{1,i};
    end
end
image = image(oversamplingLoc{1,2},oversamplingLoc{1,1},oversamplingLoc{1,3},:,:);

% turn images for correct output
if(turnImage)
    image = image(end:-1:1,:,end:-1:1,:,:);
    if(strcmp(measPara.dimension,'2D'))
        if(measPara.dim(4) > 1)
            image = image(:,:,end:-1:1,end:-1:1);
        else
            image = image(:,:,end:-1:1);
        end
    end
end

% (arbitrary) rotation correction already applied to MDH data
% % yaw (rotate around z axis)
% for i=1:size(image,3)
%     imageOut(:,:,i) = imrotate(image(:,:,i),round(2*angles(1)*180/pi), 'loose', interp);
% end
% % pitch (rotate around y axis)
% imageTmp = permute(imageOut,[2 3 1]); % => x-z-y
% clear 'imageOut';
% for i=1:size(imageTmp,3)
%     imageOut(:,:,i) = imrotate(imageTmp(:,:,i),round(angles(3)*180/pi), 'loose', interp);
% end
% % roll (rotate around x axis)
% imageTmp = permute(imageOut,[2 3 1]); % => z-y-x
% clear 'imageOut';
% for i=1:size(imageTmp,3)
%     imageOut(:,:,i) = imrotate(imageTmp(:,:,i),round(angles(2)*180/pi), 'loose', interp);
% end
% image = permute(imageOut,[2 3 1]); % => y-x-z

end