function x=nsst_rec1(dst,lpfilt)
% This function performs the inverse (local) nonsubsampled shearlet transform as given
% in G. Easley, D. Labate and W. Lim, "Sparse Directional Image Representations
% using the Discrete Shearlet Transform", Appl. Comput. Harmon. Anal. 25 pp.
% 25-46, (2008).
%
% Input:
%
% dst          - the nonsubsampled shearlet coefficients
%
% lpfilt       - the filter to be used for the Laplacian
%                Pyramid/ATrous decomposition using the codes
%                written by Arthur L. Cunha
%
% Output 
% 
% x         - the reconstructed image 
%
% Code contributors: Glenn R. Easley, Demetrio Labate, and Wang-Q Lim.
% Copyright 2011 by Glenn R. Easley. All Rights Reserved.
%

level=length(dst)-1;
y{1}=dst{1};
for i=1:level,
      y{i+1} = real(sum(dst{i+1},3));
end

x=real(atrousrec(y,lpfilt));

