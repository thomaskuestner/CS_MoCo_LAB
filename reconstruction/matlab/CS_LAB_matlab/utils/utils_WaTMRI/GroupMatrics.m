function [G,GmatZ,GmatX,groups]=GroupMatrics(parentIndex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tree structure in wavelet coefficients (2D).
% input: parentIndex. ith row are the parents of i.
%        gamma
% output: 
%         G      -- with GX=Z, X is the wavelet coefficients and Z is the
%                   coefficients with nonoverlap groups. 
%         GmatZ  -- sqrt(GmatZ*(Z^2)) is the L21 norm of Z, with group size
%                     2
%         GmatX  --  sqrt(GmatX*(X^2)) is the L21 norm of X, with group
%                     size 2. X is the original coefficients vector.
%         groups --  groups index in Z.  Z(i) belongs to the groups(i)th group.
%Author: xxxx, xxxx, 03.09.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m n]= size(parentIndex);
groupsize = n+1; 

groups = [];
for i=1:groupsize
    groups = [groups;1:m];
end
groups = groups(:);   
% 

I = 1:m;
Iy = [I' parentIndex]; %add self
Iy = Iy';
Iy = Iy(:);  % construct Y coordinates of G

zerorows = find(Iy==0);

if ~isempty(zerorows)
Iy(zerorows)=[]; %delete the parents of 0
groups(zerorows)=[]; %delete the parents of 0 
end

rowsG = size(Iy,1);
Ix = (1:rowsG)'; % X coordinates of G
G=[Ix Iy ones(rowsG,1)];
G = spconvert(G);

J=1:rowsG; 
GmatZ = [groups J'  ones(rowsG,1)];
GmatZ = spconvert(GmatZ);
GmatX = [groups Iy ones(rowsG,1)]; % P is the parent-children index in X
GmatX = spconvert(GmatX);

end