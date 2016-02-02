function [AtA,A] = corrMatrix(obj, i)

if(nargin == 2)
    [~,~,nCoil] = size(obj.kCalib{i});

    A = [];

    for n=1:nCoil
        tmp  = im2col(obj.kCalib{i}(:,:,n),obj.kernelSize,'sliding').';
        A = [A, tmp];
    end

else
    [~,~,nCoil] = size(obj.kCalib);

    A = [];

    for n=1:nCoil
        tmp  = im2col(obj.kCalib(:,:,n),obj.kernelSize(1:2),'sliding').';
        A = [A, tmp];
    end
end

AtA = A'*A;

