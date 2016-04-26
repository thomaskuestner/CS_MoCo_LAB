function [x1n,y1n,x2n,y2n,D]=gen_x_y_cordinates(n)
%
% This function generates the x and y vectors that contain
% the i,j coordinates to extract radial slices
%
% Input: n is the order of the block to be used
%
% Outputs: x1,y1 are the i,j values that correspond to 
%             the radial slices from the endpoints 1,1 to n,1
%             through the origin 
%
%           x2,y2 are the i,j values that correspond to 
%             the radial slices from the endpoints 1,1 to 1,n
%             through the origin 
%        
%           D is the matrix that contains the number of times
%           the polar grid points go through the rectangular grid 
%
% Written by Glenn Easley on December 11, 2001.
% Copyright 2011 by Glenn R. Easley. All Rights Reserved.
%
n=n+1;
% initialize vectors
x1=zeros(n,n);
y1=zeros(n,n);
x2=zeros(n,n);
y2=zeros(n,n);

xt=zeros(1,n);
m=zeros(1,n);

for i=1:n,
   y0=1; x0=i; x_n=n-i+1; y_n=n;
   if (x_n==x0),
       flag=1;
   else
      m1(i)= (y_n-y0)/(x_n-x0);
      flag=0;
   end
   xt(i,:)=linspace(x0,x_n,n);
   for j=1:n,
        if flag==0,
          y1(i,j)=m1(i)*(xt(i,j)-x0)+y0;  
          y1(i,j)=round(y1(i,j));
          x1(i,j)=round(xt(i,j));        
          x2(i,j)=y1(i,j);
          y2(i,j)=x1(i,j);
        else
          x1(i,j)=(n-1)/2+1;
          y1(i,j)=j;
          x2(i,j)=j;
          y2(i,j)=(n-1)/2+1;
        end
   end
end

n=n-1;
x1n=zeros(n,n);
y1n=zeros(n,n);
x2n=zeros(n,n);
y2n=zeros(n,n);
for i=1:n,
    for j=1:n,
         x1n(i,j)=x1(i,j);
         y1n(i,j)=y1(i,j);
         x2n(i,j)=x2(i+1,j);
         y2n(i,j)=y2(i+1,j);
     end
 end

% correct for portion outside boundry

 
x1n=flipud(x1n);
y2n(n,1)=n;
%y2n=flipud(y2n);



 D=avg_pol(n,x1n,y1n,x2n,y2n);
 
 
% end of gen_x_y_cord function

function D=avg_pol(L,x1,y1,x2,y2)
%
% This function generates the matrix that contains the number
% of times the polar grid points go through the rectangular grid 
% point i,j
%
% Input: L is the order of the block matrix
%
% Output: D is the common grid point values
%
% Written by Glenn Easley on December 11, 2001.
%


D=zeros(L);
for i=1:L,
   for j=1:L,
      D(y1(i,j),x1(i,j))=D(y1(i,j),x1(i,j))+1;
   end
end

for i=1:L,
   for j=1:L,
      D(y2(i,j),x2(i,j))=D(y2(i,j),x2(i,j))+1;
   end
end

% end of avg_pol function