function dst_scalars=nsst_scalars_e(L,shear_f,lpfilt)
% Since the nonsubsampled descrite shearlet transform is not orthogonal
% this function computes the noise level scalars of the transform with
% assigned parameters for the efficient version. 
%
% Inputs:
% 
% L                      - size of image decomposition
%
% shear_f                - the cell array containing the shearing filters 
%
%
% lpfilt                 - lpfilt is the filter to be used for the Laplacian
%                          Pyramid/ATrous decomposition using the codes
%                          written by Arthur Cunha
%
% Output:
%
% dst_scalars            - the cell array containing the scalars of 
%                          estimated noise levels for a white Gaussian noise
%                          of standard deviation 1 transform coefficients
%                          using Monte Carlo method with one iteration. 
%
% Code contributors: Glenn R. Easley, Demetrio Labate, and Wang-Q Lim.
% Copyright 2011 by Glenn R. Easley. All Rights Reserved.
%

noise=randn(L,L);
level=length(shear_f);

% LP decomposition
y_noise = atrousdec(noise,lpfilt,level);

dst_scalars=cell(1,level+1);
dst_noise=y_noise{1}; 
dst_scalars{1}=median(abs(dst_noise(:) - median(dst_noise(:))))/.6745;


for i=1:level, 
    l=size(shear_f{i},3);
    for k=1:l,          
        dst_noise=conv2(y_noise{i+1},shear_f{i}(:,:,k),'same');
        dst_scalars{i+1}(k)=median(abs(dst_noise(:) - median(dst_noise(:))))/.6745; 
    end
end
