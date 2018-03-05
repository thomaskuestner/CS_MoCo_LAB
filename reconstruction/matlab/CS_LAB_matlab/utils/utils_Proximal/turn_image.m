function [ image ] = turn_image( image )
% turn images for correct output
% code from postproc_main.m -> cropPost.m

image = image(end:-1:1,:,end:-1:1);

