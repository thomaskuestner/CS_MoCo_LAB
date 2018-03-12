function [x_cell_thresh] = softthresh_proxA_cha(x_cell,threshold,nCha,waveS,groupnorm_index)
% softthresholding combined with proximal average of wavelet tree structre
% note: calculation is equal to x/|x| (|x| - threshold)+  == x (1 - threshold/|x|)+
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

length_x = size(x_cell{1,1},2);
maxlvl_st = waveS(1,1)*waveS(1,2)+1;
maxlvl_end = 4*waveS(1,1)*waveS(1,2);
no_children = length_x-3*waveS(size(waveS,1)-1,1)*waveS(size(waveS,1)-1,2)+1;
x_cell_thresh = cell(1,nCha);
x_groupnorm_children = zeros(1,3*waveS(size(waveS,1)-1,1)*waveS(size(waveS,1)-1,2));

for j=1:nCha
    x_cell_thresh{1,j} = zeros(1,length_x);
end;
x_groupnorm = zeros(1,length_x);
% x_cell{2,j} contains also the maxlvl coeffs -> only values with children have to be calculated seperatly.
for j=1:nCha
    for n = 1:2
        x_groupnorm = x_groupnorm + x_cell{n,j}.^2;
    end;
    x_groupnorm_children = x_groupnorm_children + 2*(x_cell{1,j}(1,no_children:end).^2);
end;
threshold = threshold * sqrt(2*nCha);
x_groupnorm = sqrt(x_groupnorm);
x_groupnorm = max(1-threshold./x_groupnorm,0);
x_groupnorm_children = sqrt(x_groupnorm_children);
x_groupnorm_children = max(1-threshold./x_groupnorm_children,0);

% proximal average of groups:
for j=1:nCha
    x_cell_thresh{1,j}(1:maxlvl_st-1) = x_cell{1,j}(1:maxlvl_st-1).*x_groupnorm(1:maxlvl_st-1);
    
    x_cell_thresh{1,j}(maxlvl_st:maxlvl_end) = x_cell{1,j}(maxlvl_st:maxlvl_end).*...
        (x_groupnorm(maxlvl_st:maxlvl_end)+x_groupnorm(groupnorm_index(2,maxlvl_st:maxlvl_end))+...
        x_groupnorm(groupnorm_index(3,maxlvl_st:maxlvl_end))+x_groupnorm(groupnorm_index(4,maxlvl_st:maxlvl_end))+...
        x_groupnorm(groupnorm_index(5,maxlvl_st:maxlvl_end)))/5;
    
    x_cell_thresh{1,j}(maxlvl_end+1:no_children-1) = x_cell{1,j}(maxlvl_end+1:no_children-1).*...
        (x_groupnorm(maxlvl_end+1:no_children-1)+x_groupnorm(groupnorm_index(1,maxlvl_end+1:no_children-1))+...
        x_groupnorm(groupnorm_index(2,maxlvl_end+1:no_children-1))+x_groupnorm(groupnorm_index(3,maxlvl_end+1:no_children-1))+...
        x_groupnorm(groupnorm_index(4,maxlvl_end+1:no_children-1)))/5;
    
    x_cell_thresh{1,j}(no_children:end) = x_cell{1,j}(no_children:end).*(x_groupnorm(no_children:end)+4*x_groupnorm_children(1:end))/5;
end;
end
