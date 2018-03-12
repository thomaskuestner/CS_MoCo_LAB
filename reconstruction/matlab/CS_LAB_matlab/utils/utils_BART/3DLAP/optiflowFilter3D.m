function [uest,coeffs]=optiflowFilter3D(I1,I2,basis, K1)
%
% code modified/optimised:
% - order of indexing in the linear system matrix
%   (less shiftdim necessary)
% - matrix initialization in single precision (except filter kernel)
% 26.08.2015, Verena Neumann/Thomas Kuestner

[s1,s2,s3]=size(I1);

if length(size(basis))<=2
    [M,N]=size(basis);
    K0=round(nthroot(M,3));
    if M~=K0*K0*K0
        error('Could not identify the size of the basis filters.')
    end
    basis=reshape(basis,[K0,K0,K0,N]);
end
[K0,L0,P0,N]=size(basis);
K=(K0-1)/2;
L=(L0-1)/2;
P=(P0-1)/2;
if round(K)~=K||round(L)~=L||round(P)~=P
    error('Basis filters are not centered.')
end

calculation=2;  	% formula used to calculate velocity from filter
if nargin == 3,
    K1=2*1*K+1;         % block size first dimension
end
K2=K1;              % block size second dimension
K3=K1;              % block size third dimension

% display(['Velocity estimation with ' num2str(N) ' different ' num2str(2*K+1) 'x' num2str(2*K+1) ' filters'])
% display(['Block size = ' num2str(K1) 'x' num2str(K2) ' pixels'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% filtering with the basis %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II_1=[];
II_2=[];

% first approach -  imfilter for multi-dimensional directly
% for n=1:N
%     B1=basis(:,:,:,n);
%     B2=B1(end:-1:1,end:-1:1,end:-1:1);
%     tic;
%     J=imfilter(I2,B2,'symmetric')-imfilter(I1,B1,'symmetric');
%     disp(['Time for initial imfiltering (',num2str(n),') = ', num2str(toc)]);
%     II_1=[II_1 J(:)];
% end
% J_1=II_1(:,1);

% second approach -  use searparble filters
s=(K+2)/4;
g=@(x)exp(-x.*x/2/s^2);
K0=ceil(K);
G= g(-K0:K0);
Gd= (-K0:K0).*G;
G=G/norm(G,1);
Gd=Gd/norm((-K0:K0).*Gd,1);
Gi=G(end:-1:1); % inverse, although they are same
Gdi=Gd(end:-1:1);

% filttic=tic;
for n=1:N
    switch n
        case 1
            % Instead of using a multidimensional gaussian kernel,
            %it uses the fact that a gaussian filter can be separated in 1D gaussian kernels.
            A = imfilter(I1, reshape(G,[length(G) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(G,[length(G) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(G,[length(G) 1 1]),'symmetric');
            A = shiftdim(A,1);
            B = imfilter(I2, reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = shiftdim(B,1);
%             J1=B-A;  
        case 2
            A = imfilter(I1, reshape(Gd,[length(Gd) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(G,[length(G) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(G,[length(G) 1 1]),'symmetric');
            A = shiftdim(A,1);
            B = imfilter(I2, reshape(Gdi,[length(Gdi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = shiftdim(B,1);
%             J2=B-A;  
        case 3
            A = imfilter(I1, reshape(G,[length(G) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(Gd,[length(Gd) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(G,[length(G) 1 1]),'symmetric');
            A = shiftdim(A,1);
            B = imfilter(I2, reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gdi,[length(Gdi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = shiftdim(B,1);
%             J3=B-A;  
        case 4
            A = imfilter(I1, reshape(G,[length(G) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(G,[length(G) 1 1]),'symmetric');
            A = imfilter(shiftdim(A,1), reshape(Gd,[length(Gd) 1 1]),'symmetric');
            A = shiftdim(A,1);
            B = imfilter(I2, reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gi,[length(Gi) 1 1]),'symmetric');
            B = imfilter(shiftdim(B,1), reshape(Gdi,[length(Gdi) 1 1]),'symmetric');
            B = shiftdim(B,1);
%             J4=B-A;  
    end
    J=B-A;
    II_2=[II_2 J(:)];
end
% II_2=[J1(:) J2(:) J3(:) J4(:)];
J_2=II_2(:,1);
% disp(['Time for initial imfiltering = ', num2str(toc(filttic))]);

% evaluate their difference
% sum(sum(abs(II_1-II_2)));

% II = II_1; J = J_1; % first approach
II = II_2; J = J_2; % second approach

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% matrices needed %%%%%%%%%
%%%%%% in the linear system %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=zeros(s1*s2*s3, N-1,N-1,'single');            %single/double
b=zeros(s1*s2*s3,N-1,'single');                 %single/double
% tic
for k=1:N-1
    for l=k:N-1
        A(:,k,l) = average(II(:,k+1).*II(:,l+1));
        A(:,l,k) = A(:,k,l);
    end
    b(:,k)= -average(II(:,k+1).*J);
end
% disp(['Time for constructing A and b matrices = ', num2str(toc)]);

coeffs=zeros(s1*s2*s3,N-1,'single');            %single/double
% for k=1:s1*s2*s3
%     coeffs(k,:)=(A(:,:,k)\b(:,k))';
% end

% tic;
for k=1:(N-1)
    for l=(k+1):(N-1)
        c=A(:,l,k)./A(:,k,k);
        for m=(k+1):(N-1)
            A(:,l,m)=A(:,l,m)-c.*A(:,k,m);
        end
        A(:,l,k)=0;
        b(:,l)=b(:,l)-c.*b(:,k);
    end
end
for k=(N-1):-1:1
    coeffs(:,k)=b(:,k);
    for m=(k+1):(N-1)
        coeffs(:,k)=coeffs(:,k)-A(:,k,m).*coeffs(:,m);
    end
    coeffs(:,k)=coeffs(:,k)./A(:,k,k);
end

% disp(['Time for inversion = ', num2str(toc)]);

coeffs=[ones(s1*s2*s3,1) coeffs];

% % Exclusion of the boundaries
% p=bndindex(I1,[K+(K1-1)/2,K+(K2-1)/2, K+(K3-1)/2]);
% coeffs(p,:)=NaN;

% if nargout>2
%     err=zeros(s1,s2);
%     err(:)=sum(coeffs.*II,2);
% end

k0=(-K:K)'*ones(1,2*K+1,'single');                  %single/double
l0=k0';
k0 = repmat(k0,[1,1,(2*K+1)]);
l0 = repmat(l0,[1,1,(2*K+1)]);
p0 = shiftdim(k0,1);
basis=reshape(basis,[(2*K+1)^3 N]);
% tic;
switch calculation
    case 1
        z1=zeros(s1,s2);
        z2=zeros(s1,s2);
%         H1(K+1+(2*K+1)*K)=1;
        for n=1:N
            z1(:)=z1(:)+(basis(:,n)'*exp(-2*1i*pi*k0(:)/(2*K+1)))*coeffs(:,n);
            z2(:)=z2(:)+(basis(:,n)'*exp(-2*1i*pi*l0(:)/(2*K+1)))*coeffs(:,n);
        end
        uest=(2*K+1)/2/pi*(angle(z1./conj(z1))+1i*angle(z2./conj(z2)));
        
    case 2
        u1=zeros(s1,s2,s3,'single');                %single/double
        u11=zeros(s1,s2,s3,'single');               %single/double
        u2=zeros(s1,s2,s3,'single');                %single/double
        u22=zeros(s1,s2,s3,'single');               %single/double
        u3=zeros(s1,s2,s3,'single');                %single/double
        u33=zeros(s1,s2,s3,'single');               %single/double
        for n=1:N
            u1(:)=u1(:)-(basis(:,n)'*k0(:))*coeffs(:,n);
            u11(:)=u11(:)+sum(basis(:,n))*coeffs(:,n);
            u2(:)=u2(:)-(basis(:,n)'*l0(:))*coeffs(:,n);
            u22(:)=u22(:)+sum(basis(:,n))*coeffs(:,n);
            u3(:)=u3(:)-(basis(:,n)'*p0(:))*coeffs(:,n);
            u33(:)=u33(:)+sum(basis(:,n))*coeffs(:,n);
        end
        u_hold = [2*(u1(:)./(u11(:))), 2*(u2(:)./(u22(:))), 2*(u3(:)./(u33(:)))];
        u_hold(sqrt(sum(abs(u_hold).^2,2))>K,:) = NaN;
        uest{1} = reshape(u_hold(:,1), s1, s2, s3);
        uest{2} = reshape(u_hold(:,2), s1, s2, s3);
        uest{3} = reshape(u_hold(:,3), s1, s2, s3);
    case 3
        theta=2*pi/10;
        z1=zeros(s1,s2);
        z2=zeros(s1,s2);
%         H1(K+1+(2*K+1)*K)=1;
        for n=1:N
            z1(:)=z1(:)+(basis(:,n)'*exp(-1i*theta*k0(:)))*coeffs(:,n);
            z2(:)=z2(:)+(basis(:,n)'*exp(-1i*theta*l0(:)))*coeffs(:,n);
        end
        uest=1/theta*(angle(z1./conj(z1))+1i*angle(z2./conj(z2)));
end
% disp(['Time for extraction = ', num2str(toc)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=average(I)
K1=evalin('caller','K1');
K2=evalin('caller','K2');
K3=evalin('caller','K3');
s1=evalin('caller','s1');
s2=evalin('caller','s2');
s3=evalin('caller','s3');

% first approach
% J_1=imfilter(reshape(I,s1,s2,s3),ones(K1,K2,K3),'symmetric')/(K1*K2*K3);

II=reshape(I,s1,s2,s3);
K0=min([K1,K2,K3]);
F=ones([1 K0]);
FA = reshape(F,[length(F) 1 1]);
J_2 = imfilter(II, FA,'symmetric');
J_2 = imfilter(shiftdim(J_2,1), FA,'symmetric');
J_2 = imfilter(shiftdim(J_2,1), FA,'symmetric');
J_2 = shiftdim(J_2,1)/(K1*K2*K3);

% difference
% sum(sum(sum(abs(J_1-J_2))))

% J=J_1(:);
J=J_2(:); % second approach
return

% function p1=bndindex(I,K)
% % if K > 4,
% %     K = 4;
% % end
% 
% if length(K)==1
%     K=[K,K,K];
% end
% [M,N,P]=size(I);
% m = repmat((1:M)',[1,N,P]);
% n = repmat((1:N), [M,1,P]);
% p = repmat(shiftdim((1:P),-1), [M,N,1]);
% p1=find(~(m>=K(1)+1&m<=M-K(1)&n>=K(2)+1&n<=N-K(2)&p>=K(3)+1&p<=P-K(3)));
% return

% function [a,b]=affine(I1,I2)
% K=evalin('caller','K');
% K1=evalin('caller','K1');
% K2=evalin('caller','K2');
% s1=evalin('caller','s1');
% s2=evalin('caller','s2');
% u1=average(I1);
% u2=average(I2);
% u11=average(I1.*I1);
% u12=average(I1.*I2);
% u22=average(I2.*I2);
% d=u22-u2.*u2;
% 
% a=(u12-u1.*u2)./d;
% b=(-u12.*u2+u22.*u1)./d;
% return