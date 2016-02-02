function [AtA,A] = corrMatrix3D(obj)
% calucate 3D correlation matrix
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

nCha = size(obj.kCalib,4);


if(isreal(obj.kCalib))
    A = zeros(prod(obj.calibSize - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision);
else
    A = complex(zeros(prod(obj.calibSize - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision),zeros(prod(obj.calibSize - obj.kernelSize + 1), prod(obj.kernelSize)*nCha,obj.measPara.precision));
end
% A = [];
counter = 1;
for n=1:nCha
    if(isreal(obj.kCalib(:,:,:,n)))
        if(strcmp(obj.measPara.precision,'single'))
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colRSingle(obj.kCalib(:,:,:,n),obj.kernelSize).';
        else
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colR(obj.kCalib(:,:,:,n),obj.kernelSize).'; % before: tmp =
        end
    else
        if(strcmp(obj.measPara.precision,'single'))
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colCSingle(obj.kCalib(:,:,:,n),obj.kernelSize).';
        else
            A(:,counter:counter+prod(obj.kernelSize)-1) = im3colC(obj.kCalib(:,:,:,n),obj.kernelSize).';
        end
    end
    counter = counter + prod(obj.kernelSize);
% 	A = [A, tmp];
end

AtA = A'*A;

end