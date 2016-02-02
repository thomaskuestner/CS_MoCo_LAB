function [IdxParent, IdxChildren, Ms]=WaveRelation2D(C, S)
% WaveRelation2D: parent/children relationships for coefficients of 2-dimensional wavelet transfom
% USAGE: [IdxParent, IdxChildren, Ms]=WaveRelation2D(C, S)
% INPUT:  C: 1 x M, scaling and wavelet coefficients after 2-d wavlet transform
%         S: (N+2) x 2, 2-d wavelet decomposition structure (N denotes total decomposition levels)
%         [See instruction of WAVEDEC2 function for definition of C and S]
% OUTPUT: IdxParent: 1 x M, parent index. IdxParent(i)=j means the parent of the ith coefficient is coefficient j 
%         IdxChildren: nChildren x (M-Ms(N)), children index. IdxChildren(i,:)=[j(1),...,j(nChildren)] means
%                      the children of the ith coefficient are coefficients j(1),...,j(nChildren). Since coefficients
%                      at level N (finest level) do not have any children, they are not stored.  
%         [The order of coefficients in IdxParent and IdxChildren are the same as the coefficient order in C]
%         Ms: 1 x N, number of wavelet coefficients for each decomposition level (scaling coeff. are excluded) 

%-------------------------------
% Lihan He, ECE, Duke University
% Last change: Nov. 25, 2008
%-------------------------------

Ms=3*prod(S(2:end-1,:),2)';

Mscaling=S(1,1)*S(1,2);
% IdxParent: scaling coefficients included
IdxParent=WaveTreeParent2D(C, S);
IdxParent=[zeros(1,Mscaling),IdxParent];
idx=find(IdxParent~=0);
IdxParent(idx)=IdxParent(idx)+Mscaling;
% IdxChildren: scaling coefficients included
IdxChildren=WaveTreeChildren2D(C, S);
IdxChildren=IdxChildren+Mscaling;
IdxChildren=[zeros(4,Mscaling),IdxChildren];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IdxParent=WaveTreeParent2D(C, S)

% Find the parent index for the 2D wavelet-decomposed tree dependence.
% Consider only the wavelet coefficients, excluding the scaling coefficients.
% see instruction of MATLAB wavedec2 function for definition of C and S.
% IdxParent has the size of wavelet coefficients, which is 3*sum(prod(S(2:end-1,:),2)), 
% with IdxParent(i)=j means the parent of the ith wavelet coefficient is
% the jth wavelet coefficient.

% Decomposition level
S=S(2:end-1,:);
N=size(S,1);

% generate indicator matrix for level 2~N
% eg. if S=[8,8;16,16;32;32]  then
% indi{2}=[ 1 1 9 9 ...57 57
%           1 1 9 9 ...57 57
%           2 2 10 10...58 58
%           2 2 10 10...58 58
%           .................
%           8 8 16 16...64 64
%           8 8 16 16...64 64]
for k=2:N
    temp=reshape([1:prod(S(k-1,:))],S(k-1,:));
    temp1=zeros(S(k,1),S(k-1,2));
    temp1(1:2:end,:)=temp;
    temp1(2:2:end,:)=temp;
    indi{k}=zeros(S(k,:));
    indi{k}(:,1:2:end)=temp1;
    indi{k}(:,2:2:end)=temp1;
end

% initialization
pt_pa=0;    % parent pointer
pt_ch=0;    % children pointer
IdxParent=zeros(1,3*sum(prod(S,2)));

% find parent index
Nch=3*S(1,1)*S(1,2);  % 1 level, no parent
pt_ch=pt_ch+Nch;
%
for k=2:N
    Nch=S(k,1)*S(k,2);
    Npa=S(k-1,1)*S(k-1,2);
    %H band
    IdxParent(pt_ch+[1:Nch])=pt_pa+indi{k}(:)';
    pt_ch=pt_ch+Nch;
    pt_pa=pt_pa+Npa;
    %V band
    IdxParent(pt_ch+[1:Nch])=pt_pa+indi{k}(:)';
    pt_ch=pt_ch+Nch;
    pt_pa=pt_pa+Npa;
    %D band
    IdxParent(pt_ch+[1:Nch])=pt_pa+indi{k}(:)';
    pt_ch=pt_ch+Nch;
    pt_pa=pt_pa+Npa;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IdxChildren=WaveTreeChildren2D(C, S)

% Find the children index for the 2D wavelet-decomposed tree dependence.
% Consider only the wavelet coefficients, excluding the scaling coefficients.
% See instruction of MATLAB wavedec2 function for definition of C and S.
% IdxChildren has the size of 4 x 3*sum(prod(S(2:end-2,:),2)), where 
% 3*sum(prod(S(2:end-2,:),2)) is the number of wavelet coefficients which 
% have childrens (the finest coefficients do not have children), and each 
% parent has four children.
% IdxChildren(:,i)=[j1;j2;j3;j4] means the ith wavelet coefficient has four
% childrens, which are coefficient j1,j2,j3 and j4.

% decomposition level
S=S(2:end-1,:);
N=size(S,1);

% generate indicator matrix for level 1~N-1
% eg. if S=[8,8;16,16;32;32]  then
% indi{1}=[ 1  3  ... 237 239
%           2  4  ... 238 240
%           17 19 ... 253 255 
%           18 20 ... 254 256]
for k=1:N-1
    temp=reshape([1:prod(S(k+1,:))],S(k+1,:));
    indi{k}=zeros(4,prod(S(k,:)));
    indi{k}(1,:)=reshape(temp(1:2:end,1:2:end),[1,prod(S(k,:))]);
    indi{k}(2,:)=reshape(temp(2:2:end,1:2:end),[1,prod(S(k,:))]);
    indi{k}(3,:)=reshape(temp(1:2:end,2:2:end),[1,prod(S(k,:))]);
    indi{k}(4,:)=reshape(temp(2:2:end,2:2:end),[1,prod(S(k,:))]);
end

% initialization
pt_pa=0;    % parent pointer
pt_ch=0;    % children pointer
IdxChildren=zeros(4, 3*sum(prod(S(1:N-1,:),2)));

% find children index
Nch=3*S(1,1)*S(1,2);  % 1 level, no parent
pt_ch=pt_ch+Nch;
%
for k=1:N-1
    Npa=S(k,1)*S(k,2);
    Nch=S(k+1,1)*S(k+1,2);
    % H band
    IdxChildren(:,pt_pa+[1:Npa])=pt_ch+indi{k};
    pt_pa=pt_pa+Npa;
    pt_ch=pt_ch+Nch;
    % V band
    IdxChildren(:,pt_pa+[1:Npa])=pt_ch+indi{k};
    pt_pa=pt_pa+Npa;
    pt_ch=pt_ch+Nch;
    % D band
    IdxChildren(:,pt_pa+[1:Npa])=pt_ch+indi{k};
    pt_pa=pt_pa+Npa;
    pt_ch=pt_ch+Nch;
end

