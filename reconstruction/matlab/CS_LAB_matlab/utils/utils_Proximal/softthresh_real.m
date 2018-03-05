function x_thresh = softthresh_real(x,threshold)
% softthresholding for real valued x
% x_out = sign(x)(|x|-threshold)+ for x in R^n
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

	x_sign = sign(x);
    res = abs(x) - threshold;
    res(threshold == inf) = 0; % not necessarily needed anymore
	x_thresh = x_sign.*(res + abs(res))/2;