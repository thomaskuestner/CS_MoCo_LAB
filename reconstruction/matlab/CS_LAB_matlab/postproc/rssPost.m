function image = rssPost(imageCha,dim)
%RSS root-sum-of-squares reconstruction
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

image = cellfun(@(x) abs(x).^2, imageCha, 'UniformOutput', false);
% dim = unique(cellfun(@(x) length(size(x)), image));
image = sqrt(sum(cell2mat(shiftdim(image,-dim+1)),dim+1));


end

