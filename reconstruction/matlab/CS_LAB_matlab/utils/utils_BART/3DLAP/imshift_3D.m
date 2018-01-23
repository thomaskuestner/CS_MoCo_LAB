function J=imshift_3D(I,u,kernel)
% USAGE    : J=imshift_3D(I,u);
% FUNCTION : Shifts the 3D image I by u{1} along the lines and by u{2}
% along the columns and u{3} in the third dimension. 'u' is an arbitrary 
% vector field structured as:
%                   u{1} = u1(x,y,z)
%                   u{2} = u2(x,y,z)
%                   u{3} = u2(x,y,z). 
% If 'u' is a single vector (i.e. u1 = a, u2 = b and u3 = c) then all the 
% pixels are shifted by the same amount. Note the structure of u is alwasys
% assumed to be a cell.
%
% This function uses cubic spline interpolation used in 'interp_3D.m' with
% half point symmetric image extensions (can be changed).
%
% DATE     : 23 November 2014 (updated 28 Jan 2015 by C Gilliam)
% AUTHOR   : Thierry Blu, mailto:thierry.blu@m4x.org

[M,N,P]=size(I);

if nargin == 2,
%     kernel= 'cubicspline';
    kernel='cubicOMOMS';
% kernel = 'bilinear';
end

if length(u{1})>1               % check to see if u is a single vector or a field
    [M0,N0,P0]=size(u{1});
    if M0==M&N0==N&P0==P
        
        % Generate x, y and z:
        x = repmat((1:M)',[1,N,P]);
        y = repmat((1:N), [M,1,P]);
        z = repmat(shiftdim((1:P),-1), [M,N,1]);
        
        % Interpolation:
        J=interp_3D(double(x-u{1}),double(y-u{2}),double(z-u{3}),double(I),kernel);  %single/double
        J=cast(J,class(I));
    else
        error('Input image and flow dimensions do not match!')
    end
else
    u=[u{:}];
    integer_shift=(norm(round(u)-u)<=1e-6);

    if integer_shift
        u=round(u);
        xshift=abs(u(1));
        yshift=abs(u(2));
        zshift=abs(u(3));
        switch sign(u(1))
            case -1
                I1=[I;I(end-(1:xshift)+1,:,:)];
            case 0
                I1=I;
            case 1
                I1=[I((xshift:-1:1),:,:);I];
        end
        switch sign(u(2))
            case -1
                I1=[I1 I1(:,end-(1:yshift)+1,:)];
            case 0
            case 1
                I1=[I1(:,(yshift:-1:1),:) I1];
        end
        
        switch sign(u(3))
            case -1
                I1=cat(3,I1,I1(:,:,end-(1:zshift)+1));
            case 0
            case 1
                I1=cat(3,I1(:,:,(zshift:-1:1)), I1);
        end

        J=I1((1:M)+max(0,-u(1)),(1:N)+max(0,-u(2)), (1:P)+max(0,-u(3)));
    else
        % Generate x, y and z:
        x = repmat((1:M)',[1,N,P]);
        y = repmat((1:N), [M,1,P]);
        z = repmat(shiftdim((1:P),-1), [M,N,1]);
        
        % Interpolation:
        J=interp_3D(x-u(1),y-u(2),z-u(3),double(I),kernel);
        J=cast(J,class(I));
    end
end