function dstn=nsst_HT(dst,sigma,Tscalars,dst_scalars)
% This function performs hard thresholding of the discrete shearlet
% tranform.
%
% Input:
%
% dst          - the nonsubsampled shearlet coefficients
%
% sigma        - the standard deviation of the noise
%
% Tscalars     - a 1x3 vector containing the threshold scalars 
%                Tscalars(1) is the threshold scalar for the low-pass
%                coefficients
%                Tscalars(2) is the threshold scalar for the band-pass 
%                coefficients
%                Tscalars(3) is the threshold scalar for the high-pass
%                coefficients
%
% dst_scalars  - the cell array containing the scalars of 
%                estimated noise levels for a white Gaussian noise
%                of standard deviation 1 transform coefficients

%
% Output  
% 
% dstn         - the thresheld nonsubsampled shearlet coefficients
%
% Code contributors: Glenn R. Easley, Demetrio Labate, and Wang-Q Lim.
% Copyright 2011 by Glenn R. Easley. All Rights Reserved.
%

level=length(dst)-1;
dstn = cell(1,level+1);
dstn{1}=dst{1}.*(abs(dst{1}) > Tscalars(1)*sigma*dst_scalars{1});

for i=1:level-1,
      l=size(dst{i+1},3);
      for k=1:l, 
           dstn{i+1}(:,:,k)=dst{i+1}(:,:,k).*(abs(dst{i+1}(:,:,k)) > Tscalars(2)*sigma*dst_scalars{i+1}(k)); 
      end
end

i=level;
l=size(dst{i+1},3);
for k=1:l,
      dstn{i+1}(:,:,k)=dst{i+1}(:,:,k).*(abs(dst{i+1}(:,:,k)) > Tscalars(3)*sigma*dst_scalars{i+1}(k)); 
end