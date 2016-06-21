function [w_s]=shearing_filters_Myer(n1,level)
% This function computes the directional/shearing filters using the Meyer window
% function.
% 
% Inputs: n1 indicates the supports size of the directional filter is n1xn1
%         level indicates that the number of directions is 2^level 
%         
% Output: a sequence of 2D directional/shearing filters w_s where the
%         third index determines which directional filter is used
%
% Comments provided by Glenn R. Easley and Demetrio Labate. 
% Originally written by Glenn R. Easley on Feb 2, 2006.
% Copyright 2011 by Glenn R. Easley. All Rights Reserved.



% generate indexing coordinate for Pseudo-Polar Grid
[x11,y11,x12,y12,F1]=gen_x_y_cordinates(n1);


wf=windowing(ones(2*n1,1),2^level);
w_s=zeros(n1,n1,2^level); %initialize window array
for k=1:2^level,
    temp=wf(:,k)*ones(n1,1)';
    w_s(:,:,k)=rec_from_pol(temp,n1,x11,y11,x12,y12,F1); % convert window array into Cartesian coord.
    w_s(:,:,k)=real(fftshift(ifft2(fftshift(w_s(:,:,k)))))./sqrt(n1); 
end


