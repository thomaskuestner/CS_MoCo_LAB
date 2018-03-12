function [x_sum] = group_huber(x,threshold,nCha)
% calculates moreau-envelope for smoothed l1,2
% Marc Fischer, Feb 2016

x_group = zeros(size(x{1,1}));
for j=1:2*nCha 
    x_group = x_group + x{1,j}.^2;
end;
x_group = sqrt(x_group);
for i=1:size(x_group,2)
    if (x_group(1,i) <= threshold)
        x_group(1,i) = x_group(1,i)/(2*threshold);
    else
        x_group(1,i) = x_group(1,i)-threshold/2;
    end;
end;
x_sum = sum(x_group);
