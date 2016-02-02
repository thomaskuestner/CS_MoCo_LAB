function [x_cell_thresh] = softthresh_group_2vec(x_cell,threshold,nCha)
% softthresholding for groups (via cells)
% x_out = x/norm(x)*((norm(x)- threshold)+ for x in R^n
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

x_cell_thresh = cell(2,1);
x_groupnorm = cell(nCha,1);
for j=1:nCha
    x_groupnorm{j} = zeros(size(x_cell{1,1}));
    for n = 1:2
        x_groupnorm{j} = x_groupnorm{j} + x_cell{n,j}.^2;
    end;
    x_groupnorm{j} = sqrt(x_groupnorm{j});
    threshold = sqrt(2*nCha)*threshold;
    x_groupnorm{j} = max(1-threshold./x_groupnorm{j},0);
    
    for n = 1:2
        x_cell_thresh{n,j} = x_cell{n,j}.*x_groupnorm{j};
    end;
end;
