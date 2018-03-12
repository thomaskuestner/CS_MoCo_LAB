function [x_cell_thresh] = softthresh_group_2vec_cha(x_cell,threshold,nCha)
% softthresholding for groups (via cells)
% x_out = x/norm(x)*((norm(x)- threshold)+ for x in R^n
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

x_cell_thresh = cell(2,nCha);
x_groupnorm = zeros(size(x_cell{1,1}));
for j=1:nCha
    for n = 1:2
        x_groupnorm = x_groupnorm + x_cell{n,j}.^2;
    end;
end;
x_groupnorm = sqrt(x_groupnorm);
x_max = max(x_groupnorm-threshold,0);
x_max(threshold == inf) = 0;
x_groupnorm(x_groupnorm == 0) = 100;
x_helper = x_max./x_groupnorm;
for j=1:nCha
    for n = 1:2
        x_cell_thresh{n,j} = x_cell{n,j}.*x_helper;
    end;
end;
