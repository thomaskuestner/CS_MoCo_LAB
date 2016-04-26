function [w1] = bishrink(y1,y2,T)
% Bivariate Shrinkage Function
% Usage :
%      [w1] = bishrink(y1,y2,T)
% INPUT :
%      y1 - a noisy coefficient value
%      y2 - the corresponding parent value
%      T  - threshold value
% OUTPUT :
%      w1 - the denoised coefficient

R  = sqrt(abs(y1).^2 + abs(y2).^2);
R = R - T;
R  = R .* (R > 0);
w1 = y1 .* R./(R+T); 