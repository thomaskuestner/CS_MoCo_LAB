function [dst,shear_f]=nsst_dec1(x,shear_parameters,lpfilt)
% This function computes the (local) nonsubsampled shearlet transform as given
% in G. Easley, D. Labate and W. Lim, "Sparse Directional Image Representations
% using the Discrete Shearlet Transform", Appl. Comput. Harmon. Anal. 25 pp.
% 25-46, (2008).
%
% Inputs:
%
% x                       - input image 
%
% shear_parameters has the following fields:
%
% shear_parameters.dcomp - a vector such that .dcomp(i) indicates that the
%                          ith decomposition level has 2^decomp(i)
%                          directions. The length of the vector plus 1 is
%                          total the number of decompostions. 
%
% shear_parameters.dsize - a vector indicating the local support of the 
%                          shearing filter is .dsize(i) for 2^dcomp(i)
%                          directions. This vector is same size as .dcomp.
%
% lpfilt                 - lpfilt is the filter to be used for the Laplacian
%                          Pyramid/ATrous decomposition using the codes
%                          written by Arthur L. Cunha
%
% Output:
%
% dst                    - the cell array containing the discrete shearlet 
%                          tranform coefficients
%
% Code contributors: Glenn R. Easley, Demetrio Labate, and Wang-Q Lim.
% Copyright 2011 by Glenn R. Easley. All Rights Reserved.
%

[L,L]=size(x);
level=length(shear_parameters.dcomp);

% LP decomposition
y = atrousdec(x,lpfilt,level);

dst = cell(1,level+1);
dst{1}=y{1}; % assign low-pass coefficients to first decomposition index

shear_f=cell(1,level); % declare cell array containing shearing filters
for i=1:level, 
    w_s=shearing_filters_Myer(shear_parameters.dsize(i),shear_parameters.dcomp(i));
    for k=1:2^shear_parameters.dcomp(i),
       shear_f{i}(:,:,k) =(fft2(w_s(:,:,k),L,L)./L); 
    end
end

for i=1:level, 
    d=sum(shear_f{i},3);
    for k=1:2^shear_parameters.dcomp(i),
        shear_f{i}(:,:,k)=shear_f{i}(:,:,k)./d;
        dst{i+1}(:,:,k)=conv2p(shear_f{i}(:,:,k),y{i+1});
    end
end


