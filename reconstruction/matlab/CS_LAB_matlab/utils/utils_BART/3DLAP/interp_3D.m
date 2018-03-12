function [J,offset]=interp_3D(x,y,z,I,inttype)
% USAGE    : [J,offset]=interp(x,y,z,I,[inttype])
% FUNCTION : Interpolates the 3D image I at the (in general, noninteger) 
% positions x, y and z (usual Matlab convention for images). x, y and z are
% 3D matrices with same size. 

% DEFAULT  : The optional argument inttype can take the values
%               * 'nearestneighbor'
%               * 'bilinear' (default)
%               * 'keys'
%               * 'cubicspline'
%               * 'cubicOMOMS'
%               * 'shiftedlinear'
% Assumes images in double precision.
%
% DATE     : 23 November 2014
% AUTHOR   : Thierry Blu, mailto:thierry.blu@m4x.org
%
% modifications by Verena Neumann/Thomas Kuestner:
% as function handling takes much time, we put the cases we needed directly
% in the main function

if nargin<4
    inttype='bilinear';
end

[a,b,c]=size(I);

% Determination of the interpolation function
switch inttype
    
    case 'nearestneighbor'
        L1=-0.5;    % Support of the interpolation kernel
        L2=+0.5;
        phi=@nearest;
        
    case 'bilinear'
        L1=-1;      
        L2=+1;
        phi=@linspline;
        
    case 'keys'
        L1=-2;      
        L2=+2;
        phi=@keys;
        
    case 'cubicspline'
        L1=-2;      
        L2=+2;
        phi=@cubicspline;
        
    case 'cubicOMOMS'
        L1=-2;      
        L2=+2;
        phi=@cubicOMOMS;
        
    case 'shiftedlinear'
        tau=1/2*(1-sqrt(3)/3);
        L1=floor(-1+tau);      
        L2=ceil(1+tau);
        phi=@(x)linspline(x-tau);        
end

% Minimum and maximum row index needed in the interpolation formula
k0=floor(min(x(:))-L2+1);
k1=floor(max(x(:))-L1);
l0=floor(min(y(:))-L2+1);
l1=floor(max(y(:))-L1);
m0=floor(min(z(:))-L2+1);
m1=floor(max(z(:))-L1);

offset=[1-k0 1-l0 1-m0];

% Smallest box enclosing the image and the (x,y) positions
kk0=min(k0,1);
kk1=max(k1,a);
ll0=min(l0,1);
ll1=max(l1,b);
mm0=min(m0,1);
mm1=max(m1,c);

% Indices used in the interpolation formula
k=floor(x-L2+1);
l=floor(y-L2+1);
m=floor(z-L2+1);

%*********** Old: For 2D images **************
% % Image extension to fill the unknown pixels
% exttype='symh';                                 % options are: 
%                                                 % 'zpd' (zero-padding), 
%                                                 % 'symh' (half-point symmetry), 
%                                                 % 'symw' (whole-point symmetry), 
%                                                 % 'ppd' (periodization)
%*********************************************

% For 3D images:
% Image extension to fill the unknown pixels
exttype='symmetric';                            % options are: 
                                                % 'circular' (Pad with circular repetition of elements within the dimension), 
                                                % 'replicate' (Pad by repeating border elements of array), 
                                                % 'symmetric' (Pad array with mirror reflections of itself), 

I0=ext(I,exttype,[1-kk0 kk1-a 1-ll0 ll1-b 1-mm0 mm1-c]);
I0=I0(1-kk0+(k0:k1),1-ll0+(l0:l1), 1-mm0+(m0:m1));
[a0,b0,c0]=size(I0);
% Prefiltering when needed
switch inttype
    
    case 'cubicspline'
        J=symfilter(2/3,1/6,I);
        J=symfilter(2/3,1/6,shiftdim(J,1));
        J=symfilter(2/3,1/6,shiftdim(J,1));
        J = shiftdim(J,1);
        I0=ext(J,exttype,[1-kk0 kk1-a 1-ll0 ll1-b 1-mm0 mm1-c]);
        I0=I0(1-kk0+(k0:k1),1-ll0+(l0:l1),1-mm0+(m0:m1));

    case 'cubicOMOMS'        
        J=symfilter(13/21,4/21,I);
        J=symfilter(13/21,4/21,shiftdim(J,1));
        J=symfilter(13/21,4/21,shiftdim(J,1));
        J = shiftdim(J,1);
        I0=ext(J,exttype,[1-kk0 kk1-a 1-ll0 ll1-b 1-mm0 mm1-c]);
        I0=I0(1-kk0+(k0:k1),1-ll0+(l0:l1),1-mm0+(m0:m1));

    case 'shiftedlinear'
        z0=tau/(1-tau);
        % along 1st dimension
        I0=1/(1-tau)*filtering(1,[1 z0],I0,'causal');
        % then 2nd
        I0=1/(1-tau)*filtering(1,[1 z0],shiftdim(I0,1),'causal');
        % finally 3rd
        I0=1/(1-tau)*filtering(1,[1 z0],shiftdim(I0,1),'causal');
        % re-sort
        I0=shiftdim(I0,1);
end

% Kernel-based interpolation formula

%**** This is probably not the best approach *****
%*  The majority of the computation time is here *
J=zeros(size(x));
% tic
[DK,DL,DM]=ndgrid(0:(L2-L1-1), 0:(L2-L1-1), 0:(L2-L1-1));
for i=1:numel(DK)
    dk = DK(i); dl = DL(i); dm = DM(i);
    ind=k+dk+offset(1)+a0*(l+dl+offset(2)-1) + a0*b0*(m+dm+offset(3)-1);       % matrices are ordered along columns in Matlab
    
    % modified by VN
    % function handling takes much time! direct implementation is not as
    % clearly arranged as the one with function handles
    % but it saves about 20% of the time
    
    if strcmp(inttype,'cubicOMOMS')
        x1=1-abs(x-k-dk);
        x13=x1.*x1.*x1;
        x2=2-abs(x-k-dk);
        x23=x2.*x2.*x2;
        x2_12 = (x2>1&x2<=2);
        x2_01=(x2>0&x2<=1);
        phi1=(1/6*x23-2/3*x13+1/42*x2-2/21*x1).*x2_12+...
            (1/6*x23+1/42*x2).*x2_01;
        x1=1-abs(y-l-dl);
        x13=x1.*x1.*x1;
        x2=2-abs(y-l-dl);
        x23=x2.*x2.*x2;
        x2_12 = (x2>1&x2<=2);
        x2_01=(x2>0&x2<=1);
        phi2=(1/6*x23-2/3*x13+1/42*x2-2/21*x1).*x2_12+...
            (1/6*x23+1/42*x2).*x2_01;
        x1=1-abs(z-m-dm);
        x13=x1.*x1.*x1;
        x2=2-abs(z-m-dm);
        x23=x2.*x2.*x2;
        x2_12 = (x2>1&x2<=2);
        x2_01=(x2>0&x2<=1);
        phi3=(1/6*x23-2/3*x13+1/42*x2-2/21*x1).*x2_12+...
            (1/6*x23+1/42*x2).*x2_01;
        J=J+phi1.*phi2.*phi3.*I0(ind);
    elseif strcmp(inttype,'shiftedlinear')
        u=1-abs(x-k-dk-tau);
        phi1=max(u,zeros(size(u)));
        u=1-abs(y-l-dl-tau);
        phi2=max(u,zeros(size(u)));
        u=1-abs(z-m-dm-tau);
        phi3=max(u,zeros(size(u)));
        I0ind = I0(ind);
        J=J+phi1.*phi2.*phi3.*I0ind;
    else
%         % code for runtime analysis
%         phi1=phi(x-k-dk);
%         phi2=phi(y-l-dl);
%         phi3=phi(z-m-dm);
%         I0ind=I0(ind);
%         J=J+phi1+phi2+phi3+I0ind;

        % original code
        J=J+phi(x-k-dk).*phi(y-l-dl).*phi(z-m-dm).*I0(ind);
    end
end
% toc
% for dk=0:(L2-L1-1)
%     for dl=0:(L2-L1-1)
%         for dm=0:(L2-L1-1)
%             ind=k+dk+offset(1)+a0*(l+dl+offset(2)-1) + a0*b0*(m+dm+offset(3)-1);       % matrices are ordered along columns in Matlab
%             J=J+phi(x-k-dk).*phi(y-l-dl).*phi(z-m-dm).*I0(ind);
%         end
%     end
% end
%*************************************************
%*************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=ext(I,exttype,extsize)
[a,b,c]=size(I);
newa=a+extsize(1)+extsize(2);
newb=b+extsize(3)+extsize(4);
newc=c+extsize(5)+extsize(6);

% Extend in x:
if extsize(1)>extsize(2)
%     J=wextend(2,exttype,I,extsize(1),'bn');
    J=padarray(I,[extsize(1),0,0],exttype);
    J=J(1:newa,:,:);
else
%     J=wextend(2,exttype,I,extsize(2),'bn');
    J=padarray(I,[extsize(2),0,0],exttype);
    J=J(end+(1:newa)-newa,:,:);
end

% Extend in y:
if extsize(3)>extsize(4)
%     J=wextend(2,exttype,J,extsize(3),'nb');
    J=padarray(J,[0,extsize(3),0],exttype);
    J=J(:,1:newb,:);
else
%     J=wextend(2,exttype,J,extsize(4),'nb');
    J=padarray(J,[0,extsize(4),0],exttype);
    J=J(:,end+(1:newb)-newb,:);
end

% Extend in z:
if extsize(5)>extsize(6)
    J=padarray(J,[0,0,extsize(5)],exttype);
    J=J(:,:,1:newc);
else
    J=padarray(J,[0,0,extsize(6)],exttype);
    J=J(:,:,end+(1:newc)-newc);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=nearest(x)
u=(x<0.5&x>=-0.5);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=linspline(x)
u=1-abs(x);
max(u,zeros(size(u)));
% u=u.*(u>=0);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=keys(x)
a=-1/2;
x=abs(x);
x2=x.*x;
x3=x.*x2;
u=((a+2)*x3-(a+3)*x2+1).*(x<=1)+...
    (a*x3-5*a*x2+8*a*x-4*a).*(x>1&x<=2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=cubicspline(x)
xi=2-abs(x);
xi2=xi.*xi;
xi3=xi.*xi2;
u=(2/3-2*xi+2*xi2-1/2*xi3).*(xi>1&xi<=2)+...
    (1/6*xi3).*(xi>0&xi<=1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=cubicOMOMS(x)
x1=1-abs(x);
x13=x1.*x1.*x1;
x2=2-abs(x);
x23=x2.*x2.*x2;
u=(1/6*x23-2/3*x13+1/42*x2-2/21*x1).*(x2>1&x2<=2)+...
    (1/6*x23+1/42*x2).*(x2>0&x2<=1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=filtering(b,a,I,type)
switch type
    case 'causal'
        J=bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:,:))),I(1,:,:)*sum(b)/sum(a));
    case 'anticausal'
        I=I(end:-1:1,:,:);
        J=bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:,:))),I(1,:,:)*sum(b)/sum(a));
        J=J(end:-1:1,:,:);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=symfilter(a,b,x)
% USAGE    : y=symfilter(a,b,x)
% FUNCTION : Implements the IIR filtering with an all-pole filter with 
% z-transform of the form 1/(a+b(z+1/z)) where a and b are such that the 
% denominator has no roots on the unit circle. This is equivalent to having
% a/b complex-valued (not real) or |a/b|>2.
%
% DATE     : 19 December 2014
% AUTHOR   : Thierry Blu, mailto:thierry.blu@m4x.org

[N,K1, K2]=size(x);

% Find the root z0 that is inside the unit circle
z0=roots([b a b]);
if abs(z0(1))<1
    z0=z0(1);
else
    z0=z0(2);
end
if abs(z0)>=1
    error('The filter has poles on the unit circle!')
end
A=1/b/(z0-1/z0);

% One-pole IIR filtering of a symmetrized version of x 
z0n=filter(1,[1 -z0],[1;zeros(2*N,1)]);
y=filter(1,[1 -z0],[x;x(end:-1:1,:,:);x(1,:,:)]);
a1=(y(end,:,:)-y(1,:,:))./(1-z0n(end));

%**************************************************************************
%*****  A naive approach to N-Dim multiplication using For loops     ******
%***                Fine provided K2 is not too large                   ***
%**************************************************************************
% tic
holder = zeros(length(z0n),K1,K2);
for l1 = 1:K2,
    holder(:,:,l1) = z0n*a1(:,:,l1);
end
y=y+holder;
% toc
%**************************************************************************
%**************************************************************************


% %**************************************************************************
% %**** N-Dim multiplication using multipod matlab code by Paolo de Leva ****
% %***   Speeds up computation but can require a large amount of memory   ***
% %**************************************************************************
% % tic 
% y = y + multiprod(z0n,a1);
% % toc
% %**************************************************************************
% %**************************************************************************
% 
% 
% %**************************************************************************
% %****   N-Dim multiplication using mtimesx mex code by James Tursa     ****
% %***     Speeds up computation. Requires building mex on first use      ***
% %**************************************************************************
% % tic 
% y = y + mtimesx(z0n,a1);
% % toc
% %**************************************************************************
% %**************************************************************************

y=A*(y(1:N,:,:)+y(2*N+1-(1:N),:,:)-x);
return
