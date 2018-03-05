function [x,cost,err] = bp_ncvx(y,A,AH,lam,a,mu,Nit,pen)
% [x,cost,err] = bp_ncvx(y,A,AH,lam,a,mu,Nit,pen)
% Signal denoising with non-convex tight frame regularization
% 
% INPUT
%     y - Input signal
%     A - Tight-frame, i.e., A^H * A = rI, r > 0
%     AH - Conjugate transpose of A
%     lam - Regularization parameter
%     a - degree of non-convexity. Note: a < 1/(r*lam)
%     mu - Augmented Lagrangian parameter for ADMM (mu > 1)
%     Nit - Number of Iterations
%     pen - Regularizer to be used 
%             a. Logarithmic ('log')
%             b. Rational    ('rat')
%             c. Arctangent  ('atan')
%             d. L1 norm     ('l1')
%    
% OUTPUT           
%     x - Denoised signal
%     cost - Cost function history
%     err - Error when using variable splitting, i.e., ||u-Ax||_2^2
%     
% Ankit Parekh (ankit.parekh@nyu.edu), NYU School of Engineering
% Reference:
% Convex denoising using non-convex tight frame regularization
% Ankit Parekh and Ivan W. Selesnick. 
% IEEE Signal Process. Lett., 2015    

%Initialize
y = y(:);
N = length(y); 
x = zeros(N,1);
d = zeros(size(A(y)));
k = 1/(1 + mu);
cost = zeros(1,Nit);
err = zeros(1,Nit);
switch pen
    case 'log'
        phi = @(x, a) (1/a) * log(1 + a*abs(x));  
    case 'rat'
        phi = @(x, a) abs(x)./ (1 + (a/2) * abs(x));
    case 'atan'
        phi = @(x, a) 2./(a*sqrt(3)) .* (atan((2*a.*abs(x)+1)/sqrt(3)) - pi/6);
    case 'l1'
        phi = @(x,a) abs(x);
end

%Iterative thresholding until convergence or user defined iterations
for i = 1:Nit
    u = thresh(A(x) + d, lam/mu,a,pen);
    x = k * (y + mu * AH(u-d));
    d = d-(u-A(x));
    
    %Calculate cost function
    residual = y-x;
    Ax = A(x);
    cost(i) = 0.5*sum(abs(residual(:)).^2) + lam*sum(phi(Ax(:),a));
    e = u-Ax;
    err(i) = sum(abs(e(:)).^2);
end

