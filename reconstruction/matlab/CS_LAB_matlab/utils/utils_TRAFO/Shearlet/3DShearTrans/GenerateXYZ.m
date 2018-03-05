%% function [X,Y,Z,D]= GenerateXYZ(n)
%% Generates the XYZ coordiante of Polar Grid in Cartesian Grid
%%Input: 
%%      n: size of the cubic box .
%%Output: X,Y,Z: coordinate of polar grid in cartesian grid.
%%        i.e X(i,k,j),Y(i,k,j),Z(i,k,j) contain cartesian coordinate of
%%         polar point(i,k,j)
%%        D: number of time polar grid point go rhough same cartesian grid
%%       point
%% 

function [X,Y,Z]= GenerateXYZ(n)
X=zeros(n,n,n,'uint16');
Y=zeros(n,n,n,'uint16');
Z=zeros(n,n,n,'uint16');
midPoint=ceil(n/2);
for j=1:midPoint  
    for i=1:midPoint
      %join the point(1,i,j)   (n, n-i+1,n-j+1)
      % create n point between them
      diagonalL2=sum(([n   n-i+1 ]-[1 i]).^2);
      diagonalL3= (diagonalL2+ (n-j+1-j)^2)^.5;
      diagonalL2=diagonalL2^.5;
      sinAlpha= (n-j+1-j)/diagonalL3;
      cosAlpha=diagonalL2/diagonalL3;
      sinBeta= (n-i+1-i) /diagonalL2;
      cosBeta= (n-1)/diagonalL2;
      deltaL3=diagonalL3/(n-1);
      for k=1:n
        %create n point along diagonal and store their X,Y,Z        
        l=deltaL3*(k-1);
        temp=l*sinAlpha;
        Z(i,k,j)=j+temp;   
        Z(n-i+1,k,j)=Z(i,k,j);
        Z(i,k,n-j+1)=n-j+1-temp;
        Z(n-i+1,k,n-j+1)=n-j+1-temp;
        
        temp=l*cosAlpha;
        tmep1=temp*sinBeta;
        Y(i,k,j)=i+tmep1;
        Y(n-i+1,k,j)=n-i+1-tmep1;
        
        X(i,k,j)=1+temp*cosBeta;   
        X(n-i+1,k,j)=X(i,k,j);
        %symmtry across lines via i
        %handle symmetry across plane via j        
        X(i,k,n-j+1)=X(i,k,j);
        X(n-i+1,k,n-j+1)=X(i,k,j);
        Y(i,k,n-j+1)= Y(i,k,j);
        Y(n-i+1,k,n-j+1)= Y(n-i+1,k,j);
        
      end
    end 
end
%F=PolFreq(n,X,Y,Z);
%% function function D=Avg(n,X,Y,Z)
%% number of time polar grid point go rhough same cartesian grid
%%       point
%%Input: 
%%     n: size of the cubic box .
%%     X,Y,Z: coordinate of plolar grid
%% Ouput:
%%       D: number of time polar grid point go rhough same cartesian grid
%%       point
%% 
% function F=PolFreq(n, X, Y,Z)
% F=zeros(n,n,n,'int16');
% for j =1:n
%   for i=1:n
%     for k=1:n
%     F(X(i,k,j),Y(i,k,j),Z(i,k,j))=F(X(i,k,j),Y(i,k,j),Z(i,k,j))+1;
%     end
%   end
% end
