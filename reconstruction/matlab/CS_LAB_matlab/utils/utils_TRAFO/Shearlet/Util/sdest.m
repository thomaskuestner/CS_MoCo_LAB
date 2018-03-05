% Estimation of the standard deviation of an image assuming that the noise
% is Gaussian white additive. The file estimates the standard deviation
% using the median filter on the fine scale subband of the wavelet
% decomposition.
 function sd_estimate=sdest(x)
 [ca,ch,cv,cd] = dwt2(x,'sym4','mode','sym');
 sd_estimate = mad(cd(:),1)/0.6745


