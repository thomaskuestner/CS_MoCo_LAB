function res = im2row(im, winSize)
%res = im2row(im, winSize)
[sx,sy,sz] = size(im);

res = zeros((sx-winSize(1)+1)*(sy-winSize(2)+1),prod(winSize),sz);
count=0;
for y=1:winSize(2)
    for x=1:winSize(1)
        count = count+1;
        res(:,count,:) = reshape(im(x:sx-winSize(1)+x,y:sy-winSize(2)+y,:),...
            (sx-winSize(1)+1)*(sy-winSize(2)+1),1,sz);
    end
end
