function I=GetToyImage(F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% force the wavelet coeffieciets of input image to fit tree structure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 0  % group-wise
    fprintf(1,'Error: Please input an image!');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D wavelet transformation and sparse signal -- theta0 
DecLevel=4;
%WaveletName='db1';
[m,n]=size(F);f=F(:);N=m*n;
[IdxParent, IdxChildren, Ms] = WaveTreeStructure2D(m,DecLevel);
wav = daubcqf(2);
WT = @(x) reshape(mdwt(reshape(x,m,n),wav,DecLevel),N,1);W = @(x) reshape(midwt(reshape(x,m,n),wav,DecLevel),N,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



theta = WT(f); %wavelet coe
theta2=theta.*theta;
st = sort(theta2);
sum2=sum(theta'*theta);

threshhold = 0.001*sum2;
t = 0;
for i=1:N
    t=t+st(i);
    if t>threshhold
        break;
    end
end

threshhold2 = st(i);
zeorrows=find(theta2<threshhold2);
theta(zeorrows)=zeros(length(zeorrows),1);


% Get mean value of last 2 levels coefficients
% thetaImg = reshape(theta,m,n);
% sumc = sum(abs(theta));
% sumlevels = sum(sum(abs(thetaImg(1:64,1:64)))); % first 2 levels
% sumc = sumc - sumlevels;
% meanc = sumc/(Ms(3)+Ms(4));
%figure;imshow(reshape(W(theta),m,n),[]);
%Z=G*theta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if parent is zero, set children to be zeros.
parent0 = find(theta==0);  %parent nodes with 0 values
par = IdxChildren(:,1);  %nodes who have children
for i=1:length(parent0)
    indexparent = find(par==parent0(i));
    parentlevel = GetLevel(parent0(i),Ms,n);
    if parentlevel<3
        continue;
    end
    if length(indexparent)==1
        children = IdxChildren(indexparent,2:5);
        for j=1:4
            if children(j)~=0 %&& (abs( theta(children(j)))<meanc)
                theta(children(j))=0;
            end
        end
    end
end


% if child is zero, set parent to be zeros.
parent = IdxParent(parent0);
parent(find(parent==0))=[];
for i=1:length(parent)
    indexparent = parent(i);
    parentlevel = GetLevel(indexparent,Ms,n);
     if parentlevel<3
            continue;
     end
     %if abs(theta(parent(i))<meanc)
         theta(parent(i))=0;
    % end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l=0;u=255;
if((l==-Inf)&&(u==Inf))
    project=@(x)x;
elseif (isfinite(l)&&(u==Inf))
    project=@(x)(((l<x).*x)+(l*(x<=l)));
elseif (isfinite(u)&&(l==-Inf))
     project=@(x)(((x<u).*x)+((x>=u)*u));
elseif ((isfinite(u)&&isfinite(l))&&(l<u))
    project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
else
    error('lower and upper bound l,u should satisfy l<u');
end
% zerorows = find(par==0);
% par(zerorows)=[];
% wf2(par)=0;
 I = reshape(W(theta),m,n);
 I=project(I);
 %imwrite(uint8(I),'toyimage.jpg');
% figure; imshow(F,[]);
 
end