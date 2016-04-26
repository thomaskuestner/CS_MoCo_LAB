function [y] = symextend(x,Nnum)

y(:,1:Nnum) = fliplr(x(:,1:Nnum));
y = [fliplr(x(:,1:Nnum)) x x(:,end:-1:end-Nnum+1)];
y = [flipud(y(1:Nnum,:)); y ;y(end:-1:end-Nnum+1,:)];

