function [x_thresh_vec1, x_thresh_vec2] = groupthresh_2vec(x_vec1,x_vec2,threshold)
% softthresholding for groups (via vectors)
% x_out = x/norm(x)*((norm(x)- threshold)+ for x in R^n
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

threshold = threshold * sqrt(1/2);
x_groupnorm = x_vec1.^2 + x_vec2.^2;
x_groupnorm = sqrt(x_groupnorm);
x_max = max(x_groupnorm-threshold,0);

% make sure values stay in 0,1
x_max(threshold == inf) = 0;
x_groupnorm(x_groupnorm == 0) = 1;
x_helper = x_max./x_groupnorm;

x_thresh_vec1 = x_vec1.*x_helper;
x_thresh_vec2 = x_vec2.*x_helper;