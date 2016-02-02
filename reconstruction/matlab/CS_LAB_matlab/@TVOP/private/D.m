function res = D(image)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005

[sx,sy] = size(image);

Dx = image([2:end,end],:) - image;
Dy = image(:,[2:end,end]) - image;

%res = [sum(image(:))/sqrt(sx*sy); Dx(:);  Dy(:)]; 
res = cat(3,Dx,Dy);


