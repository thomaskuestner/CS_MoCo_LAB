function C=rec_from_pol(l,n,x1,y1,x2,y2,D)
%
% This funcion re-assembles the radial slice into block.
%
% Inputs: l is the radial slice matrix
%         n is the order of the matrix that is to be re-assembled
%         x1,y1,x2,y2 are the polar coordinates generated from function
%                gen_x_y_cord
%         D is the matrix containing common polar grid points
% Output: C is the re-assembled block matrix 
%
% Written by Glenn Easley on December 11, 2001.
% Copyright 2006 by Glenn R. Easley. All Rights Reserved.
%

C=zeros(n,n);
option=0;
if (option==1),
for i=1:n,
    for j=1:n,
       C(y1(i,j),x1(i,j))=l(i,j);
       C(y2(i,j),x2(i,j))=l(i+n,j);
    end
 end
else
for i=1:n,
    for j=1:n,
       C(y1(i,j),x1(i,j))=C(y1(i,j),x1(i,j))+l(i,j);
       C(y2(i,j),x2(i,j))=C(y2(i,j),x2(i,j))+l(i+n,j);
    end
 end
 
% average common radial grid values 
C=C./D;
end
% end of rec_from_pol function



