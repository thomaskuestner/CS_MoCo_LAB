function [IdxParent, IdxChildren, Ms] = WaveTreeStructure2D(N, L);
% WaveTreeStructure: parent/children relationships for coefficients 
%     of 2-dimensional wavelet transfom
% Invoked after W = @(x) midwt(x,wav);WT = @(x) mdwt(x,wav);
% It finds the relations between the coefficient images obtained by Rice
%     wavelet toolbox
%
% USAGE: [IdxParent, IdxChildren]=WaveTreeStructure2D(C, S)
%
% INPUT:  N: image with size of N x N, scaling/coefficients 
%            after 2-d wavelet transform "mdwt"
%         L: denotes total decomposition levels
%         [See instruction of "mdwt" function for definition of C and L]
%
% OUTPUT: IdxParent: N^2 x 1, parent index. 
%         IdxParent(i)=j means the parent of the ith coefficient is coefficient j 
%
%         IdxChildren: nChildren x (M-Ms(N)), children index. 
%         IdxChildren(i,2:end)=[j(1),...,j(nChildren)] means the children of 
%         the ith coefficient are coefficients j(1),...,j(nChildren). 
%         Since coefficients at level L (finest level) do not have any children, 
%         they are not stored.  
%
%         The order of coefficients in IdxParent and IdxChildren are the 
%         same as the coefficient order in C(:)
%
%         Ms: 1 x N, number of wavelet coefficients for each decomposition 
%             level (scaling coeff. are excluded) 

%--------------------------------------------------------------------------
% Junzhou Huang, CSE, University of Texas at Arlington, Feb, 2012
%--------------------------------------------------------------------------
for i=1:L    
    Ms(i,1)=3*(N/2^(L-i+1))^2;
end
%---------------------------------------------------------------------------%
% INITIALIZATIONS
IdxParent=zeros(N^2,1);
IdxChildren=zeros(N^2/4,5); %Col 1: index i; col 2-5 the 4 children of index i

IJ=zeros(N^2, 2);
%    col 1: R  (row pointer)
%    col 2: C  (column pointer)
[Isub, Jsub]=ind2sub([N, N], [1:N^2]');
IJ(:,1)=Isub;IJ(:,2)=Jsub;

%%%% Parent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=zeros(N, N); 
B(1:N/2^(L-1),1:N/2^(L-1)) = ones(N/2^(L-1), N/2^(L-1));
index0=find(B==1);
IdxParent(index0)=0;

index=find(B==0);
[Isub, Jsub]=ind2sub([N, N], index);

indexPar=sub2ind([N, N], ceil(Isub/2), ceil(Jsub/2));
IdxParent(index)=indexPar;

%%%%%%%% child  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=zeros(N, N); 
B(1:N/2^(L-1),1:N/2^(L-1)) = ones(N/2^(L-1), N/2^(L-1));
index0=find(B==1);

IdxChildren(1:length(index0), 1)=index0;
IdxChildren(1:length(index0), 2:5)=0;

B(N/2+1:N, :) = ones(N/2, N);
B(:, N/2+1:N) = ones(N, N/2);

index=find(B==0);
[Isub, Jsub]=ind2sub([N, N], index);
indexChild4=sub2ind([N, N], Isub*2, Jsub*2);
indexChild3=sub2ind([N, N], Isub*2, Jsub*2-1);
indexChild2=sub2ind([N, N], Isub*2-1, Jsub*2);
indexChild1=sub2ind([N, N], Isub*2-1, Jsub*2-1);

IdxChildren(length(index0)+1:end,:)=[index, indexChild1, indexChild2, indexChild3, indexChild4];
IdxChildren2=IdxChildren;
[YY,II] = sort(IdxChildren(:,1));

IdxChildren=zeros(size(IdxChildren));
IdxChildren=IdxChildren2(II,:);



















