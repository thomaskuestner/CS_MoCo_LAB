function [x,cost] = bp_ncvxUDWT(y,A,AH,J,lam,a,mu,Nit,pen)
% 
% function [x,cost] = bp_ncvxUDWT(y,A,AH,J,lam,a,mu,Nit,pen)
% 1D signal denoising using non-convex regularization with the undecimated
% wavelet transform
% Input:  
%     y - Input Signal
%     A - Forward transform (Undecimated Wavelet transform)
%     AH - Inverse transform
%     J - Number of scales
%     lam - Regularization parameter (vector)
%     a - Degree of non-convexity (a_i < 1/lam_i)
%     mu - Augmented Lagrangian parameter (mu > 1/r)
%     Nit - Number of iterations
%     pen - Regularizer to be used 
%             a. Logarithmic ('log')
%             b. Rational    ('rat')
%             c. Arctangent  ('atan')
%             d. L1          ('l1')
%     
% Output:
%     x - Denoised Signal
%     cost - Cost function history
%         
% Ankit Parekh (ankit.parekh@nyu.edu), NYU School of Engineering
% Reference:
% Convex denoising using non-convex tight frame regularization
% Ankit Parekh and Ivan W. Selesnick
% IEEE Signal Process. Lett., 2015

%Initialization
x = zeros(size(y));
d = A(y);
k = 1/(1 + mu);
cost = zeros(1,Nit);

switch pen
    case 'log'
        phi = @(x, r) (1/r) * log(1 + r*abs(x));                    
    case 'atan'
        phi = @(x, a) 2./(a*sqrt(3)) .* (atan((2*a.*abs(x)+1)/sqrt(3)) - pi/6);
    case 'rat'
        phi = @(x, a) abs(x)./ (1 + (a/2) * abs(x));
    case 'l1'
        phi = @(x,a) abs(x);
end

%Iterative procedure
for i = 1:Nit
    % Thresholding step
    w = A(x);
    for j = 1:J
        C = (w{j} + d{j});
        w{j} = thresh(C,lam{j}/mu,a{j},pen);
    end
    
    % Compute iterate    
    for j = 1:J
       U{j} = w{j} - d{j};     
    end
    U{J+1} = w{J+1}; 
    x = k * (y + mu * AH(U));
    u = A(x);
    
    % Update splitting
    for j = 1:J
        d{j} = d{j} - (w{j} - u{j});
    end
    
    %Calculate cost function
    residual = y-x;
    Ax = A(x);
    P = 0;
    for r = 1:J
        P = P + sum(lam{r}.*phi(Ax{r}, a{r}));
    end
    cost(i) = 0.5*sum(abs(residual(:)).^2) + P;
   
end

