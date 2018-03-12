function [AtA,A] = corrMatrix4D(obj,i)
% calucate 4D correlation matrix
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

nCha = size(obj.kCalib{i},4);

% A = [];

if(isreal(obj.kCalib))
    A = zeros(prod(obj.calibSize{i} - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision);
else
    A = complex(zeros(prod(obj.calibSize{i} - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision),zeros(prod(obj.calibSize{i} - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision));
end
if(isempty(A))
    error('corrMatrix4D(): Unable to create GRAPPA kernel. Check kernel and calibration size dimensionality');
end
counter = 1;
for n=1:nCha
    if(isreal(obj.kCalib{i}(:,:,:,n)))
        if(strcmp(obj.measPara.precision,'single'))
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colRSingle(obj.kCalib{i}(:,:,:,n),obj.kernelSize).';
        else
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colR(obj.kCalib{i}(:,:,:,n),obj.kernelSize).'; % before: tmp =
        end
    else
        if(strcmp(obj.measPara.precision,'single'))
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colCSingle(obj.kCalib{i}(:,:,:,n),obj.kernelSize).';
        else
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colC(obj.kCalib{i}(:,:,:,n),obj.kernelSize).';
        end
    end
    counter = counter + prod(obj.kernelSize);
end

AtA = A'*A;

end