function J=ShiftedLinear_Interp_3D(ux,uy,uz,I)
% USAGE    : J=ShiftedLinear_Interp_3D(ux,uy,uz,I)
% FUNCTION : Interpolates the 3D image I at the (in general, noninteger) 
% positions x - ux, y - uy and z - uz. ux, uy and uz are 3D matrices with 
% same size and represent displacements in x, y and z, respectively. The 
% function combines shifted-linear interpolation with matlab's inbuilt 
% interp3 (with Trilinear input). The result is a fast interpolation that 
% is more accurate than standard Trilinear interpolation (approximately 
% equivalent to Key's cubic interpolation). For more details see:
%
% T. Blu, P. Thevenaz, and M. Unser, "Linear interpolation revitalized,"
%   IEEE Trans. Image Processing, vol. 13, no. 5, pp. 710–719, 2004.
%
% DATE     : 10 March 2016
% AUTHOR   : Christopher Gilliam, email: dr.christopher.gilliam@ieee.org
%
% Note: 
% 1) Based on the more general imshift/interp function written by 
%   Thierry Blu, mailto:thierry.blu@m4x.org. 
%
% 2) The linear interpolation is achieved using the mex file written by 
%   Andriy Myronenko. It is faster than matlab's inbuilt interp3 function.
% 

[a,b,c]=size(I);

% Generate x, y and z:
x = repmat((1:a)',[1,b,c]);
y = repmat((1:b), [a,1,c]);
z = repmat(shiftdim((1:c),-1), [a,b,1]);

% Parameters for shifted-linear interpolation:
tau=1/2*(1-sqrt(3)/3);
L1=floor(-1+tau);      
L2=ceil(1+tau);

x1 = double(x - ux);
y1 = double(y - uy);
z1 = double(z - uz);

% Minimum and maximum row index needed in the interpolation formula
k0=floor(min(x1(:))-L2+1);
k1=floor(max(x1(:))-L1);
l0=floor(min(y1(:))-L2+1);
l1=floor(max(y1(:))-L1);
m0=floor(min(z1(:))-L2+1);
m1=floor(max(z1(:))-L1);

% Smallest box enclosing the image and the (x,y) positions
kk0=min(k0,1);
kk1=max(k1,a);
ll0=min(l0,1);
ll1=max(l1,b);
mm0=min(m0,1);
mm1=max(m1,c);


% For 3D images:
% Image extension to fill the unknown pixels
exttype='symmetric';                            % options are: 
                                                % 'circular' (Pad with circular repetition of elements within the dimension), 
                                                % 'replicate' (Pad by repeating border elements of array), 
                                                % 'symmetric' (Pad array with mirror reflections of itself), 

I0=ext(I,exttype,[1-kk0 kk1-a 1-ll0 ll1-b 1-mm0 mm1-c]);
I0=I0(1-kk0+(k0:k1),1-ll0+(l0:l1), 1-mm0+(m0:m1));
[a0,b0,c0]=size(I0);
x = repmat((k0:k1)',[1,b0,c0]);
y = repmat((l0:l1), [a0,1,c0]);
z = repmat(shiftdim((m0:m1),-1), [a0,b0,1]);

% Prefiltering image I0 so as to have shifted-linear interpolation:
z0=tau/(1-tau);
% along 1st dimension
I0=1/(1-tau)*filtering(1,[1 z0],I0,'causal');
% then 2nd
I0=1/(1-tau)*filtering(1,[1 z0],shiftdim(I0,1),'causal');
% finally 3rd
I0=1/(1-tau)*filtering(1,[1 z0],shiftdim(I0,1),'causal');
% re-sort
I0=shiftdim(I0,1);

% % Matlab's inbuilt interpolation command:
% J = interp3(y,x,z,I0,y1-tau,x1-tau,z1-tau,'linear');

% Andriy Myronenko's interpolation code:
y1 = y1 - (min(y(:))-1);
x1 = x1 - (min(x(:))-1);
z1 = z1 - (min(z(:))-1);
J = mirt3D_mexinterp(double(I0), y1-tau, x1-tau,z1-tau);

J=cast(J,class(I));

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