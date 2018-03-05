function F=PolFreq(n, X, Y,Z)
F=zeros(n,n,n);
for j =1:n
  for i=1:n
    for k=1:n
    F(X(i,k,j),Y(i,k,j),Z(i,k,j))=F(X(i,k,j),Y(i,k,j),Z(i,k,j))+1;
    end
  end
end