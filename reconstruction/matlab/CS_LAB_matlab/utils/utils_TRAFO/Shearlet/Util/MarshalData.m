function X= MarshalData(X,level)

divFactor=int8([3 6 12 24]);
[x, y ,z]=size(X);

xRem= divFactor(level)- rem(x,divFactor(level));
yRem= divFactor(level)- rem(y,divFactor(level));
zRem= divFactor(level)- rem(z,divFactor(level));

% Y=zeros(x+xRem,y+yRem,z+zRem);
% Y(1:x,1:y,1:z)=X(1:x,1:y,1:z);
for i=1:xRem
X(end+1,:,:)=X(x,:,:);
end

for i=1:yRem
X(:,end+1,:)=X(:,y,:);
end

for i=1:zRem
X(:,:,end+1)=X(:,:,z);
end

end