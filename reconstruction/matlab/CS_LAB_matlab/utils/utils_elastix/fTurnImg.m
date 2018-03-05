function dImg = fTurnImg( dImg, lDirection, iPerm )
%FTURNIMG turn image (4D, RSS channel combined or not) into coronal view for MC 
% (lDirection = 0) or backwards to CS recon dimension (lDirection = 1)
% ndims == 5
% lDirection = 0: t-y-z-x-cha => x-y-z-t-cha
% lDirection = 1: x-y-z-t-cha => t-y-z-x-cha
% ndims == 4
% lDirection = 0: t-y-z-x => x-y-z-t
% lDirection = 1: x-y-z-t => t-y-z-x
% ndims == 3
% lDirection = 0: NOP
% lDirection = 1: x-y-z => y-z-x 
% for all other cases define the permutation vector in iPerm
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

% mapping the input dimensions
if(nargin < 3)
    if(lDirection == 0)
        if(ndims(dImg) == 5)
            iPerm = [4, 2, 3, 1, 5]; % => x-y-z-t-cha
        elseif(ndims(dImg) == 4)
            iPerm = [4, 2, 3, 1]; % => x-y-z-t
        else
            error('fTurnImg(): Unknown image dimensionality');
        end
    elseif(lDirection == 1)
        if(ndims(dImg) == 5)
            iPerm = [4, 2, 3, 1, 5]; % => t-y-z-x-cha
        elseif(ndims(dImg) == 4)
            iPerm = [4, 2, 3, 1]; % => t-y-z-x-cha
        elseif(ndims(dImg) == 3)
            iPerm = [2, 3, 1];
        else
            error('fTurnImg(): Unknown image dimensionality');
        end
    end
end

if(lDirection == 0)
    dImg = permute(dImg, iPerm);
    dImg = flipdim(dImg, 1);
    dImg = flipdim(dImg, 2);
%     dImg = flipdim(dImg, 4);
elseif(lDirection == 1)
%     dImg = flipdim(dImg, 4);
    dImg = flipdim(dImg, 2);
    dImg = flipdim(dImg, 1); 
    dImg = permute(dImg, iPerm); 
end

