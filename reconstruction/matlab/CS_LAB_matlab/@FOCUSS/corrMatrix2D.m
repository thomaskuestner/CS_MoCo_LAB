function [AtA,A] = corrMatrix2D(obj,i)
% calucate 2D correlation matrix
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

nCha = size(obj.kCalib{i},3);

% A = [];
if(isreal(obj.kCalib{i}))
    A = zeros(prod(obj.calibSize - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision);
else
    A = complex(zeros(prod(obj.calibSize - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision),zeros(prod(obj.calibSize - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision));
end

counter = 1;
for n=1:nCha
	A(:,counter:counter+prod(obj.kernelSize)-1)  = im2col(obj.kCalib{i}(:,:,n),obj.kernelSize,'sliding').';
    counter = counter + prod(obj.kernelSize);
% 	A = [A, tmp];
end

AtA = A'*A;

end