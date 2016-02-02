function res = wavMask(imSize,scale)
% function scales the value of each wavelet scale such they display nicely.

sx = imSize(1);
sy = imSize(2);
res = zeros(imSize)+1;
NM = round((log2(imSize)));
for n = 1:min(NM)-scale+1
	res(1:round(2^(NM(1)-n)),1:round(2^(NM(2)-n))) = res(1:round(2^(NM(1)-n)),1:round(2^(NM(2)-n)))/2;
end
%res = sqrt(res);
%res = res/min(res(:));
%n = n+1;
%res(1:round(2^(NM(1)-n)),1:round(2^(NM(2)-n)))=0;

	

