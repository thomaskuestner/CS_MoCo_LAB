function res = adjD(y)

res = zeros(size(y,1),size(y,2));

%y1 = ones(imsize)*y(1)/sqrt(prod(imsize));
%yx = (reshape(y(2:prod(imsize)+1), imsize(1), imsize(2)));
%yy = (reshape(y(prod(imsize)+2:end), imsize(1), imsize(2)));

res = adjDx(y(:,:,1)) + adjDy(y(:,:,2));

return;


function res = adjDy(x)
res = x(:,[1,1:end-1]) - x;
res(:,1) = -x(:,1);
res(:,end) = x(:,end-1);

function res = adjDx(x)
res = x([1,1:end-1],:) - x;
res(1,:) = -x(1,:);
res(end,:) = x(end-1,:);


