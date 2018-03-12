function scaledImg = scaleImg(img,range)
% scale image to range
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

if(nargin < 2)
    range = [0 1];
end

if(length(range) == 1)
    range = [0 range];
end

scaledImg = ((img - min(img(:))) * (range(2)-range(1)))./(max(img(:)) - min(img(:)));

end