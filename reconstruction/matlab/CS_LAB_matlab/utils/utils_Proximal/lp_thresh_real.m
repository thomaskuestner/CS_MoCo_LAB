function [ y ] = lp_thresh_real( x,threshold,p )
% generalized softthresholding (GST) based on the paper by Zuo et al.
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

threshold_nonconvex = (2*threshold*(1-p))^(1/(2-p)) + threshold*p*(2*threshold*(1-p))^((p-1)/(2-p));
y = zeros(size(x));
i0 = find(abs(x)>threshold_nonconvex);

if p < 1
    J = 3;
else
    J = 1;
end;
% if length(i0) > 1
    x0 = x(i0);
    x_temp = abs(x0);
    for j = 1:J
        x_temp = abs(x0) - threshold*p*(x_temp.^(p-1));
    end;
    y(i0) = sign(x0).*x_temp;
% end