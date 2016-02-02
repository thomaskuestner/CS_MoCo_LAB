function [AtA,A] = corrMatrix(kCalib, kSize)


disp('corrMatrix.m Warning: This function is deprecated. Use dat2AtA instead')

[sx,sy,nCoil] = size(kCalib);

A = [];
y = [];


for n=1:nCoil
	tmp  = im2col(kCalib(:,:,n),kSize,'sliding').';
	A = [A, tmp];
end


AtA = A'*A;

