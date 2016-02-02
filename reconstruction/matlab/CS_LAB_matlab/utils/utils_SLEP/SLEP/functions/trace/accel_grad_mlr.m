% revised on September 21, 2010 to make parameters consistent

% September 25, 2009
% written by Shuiwang Ji and Jieping Ye

% This function implements the accelerated gradient algorithm for 
% multivariate linear regression regularized by trace norm as
% described in Ji and Ye (ICML 2009).

% References:
%Ji, S. and Ye, J. 2009. An accelerated gradient method for trace norm minimization. 
%In Proceedings of the 26th Annual international Conference on Machine Learning 
%(Montreal, Quebec, Canada, June 14 - 18, 2009). ICML '09, vol. 382. ACM, New York, 
%NY, 457-464.

%[W, fval_vec, itr_counter] = accel_grad_mlr(A,Y,lambda,opt)

% required inputs:
% A: N x D where each each row is a data point and N is the sample size
% Y: N x M where M is the number of regression tasks
% lambda: regularization parameter

% optional inputs:
% opt.L0: Initial guess for the Lipschitz constant
% opt.gamma: the multiplicative factor for Lipschitz constant
% opt.W_init: initial weight matrix
% opt.epsilon: precision for termination
% opt.max_itr: maximum number of iterations

% outputs:
% W: the computed weight matrix
% fval_vec: a vector for the sequence of function values
% itr_counter: number of iterations executed

function [W, fval_vec, itr_counter] = accel_grad_mlr(A,Y,lambda,opt)

Xtrain = A;
Ytrain = Y;
clear A;
clear Y;

if nargin<4
    opt = [];
end

if isfield(opt, 'L0')
    L0 = opt.L0;
else
    L0 = 100;
end

if isfield(opt, 'gamma')
    gamma = opt.gamma;
else
    gamma = 1.1;
end

if isfield(opt, 'W_init')
    W_init = opt.W_init;
else
    W_init = zeros(size(Xtrain,2),size(Ytrain,2));
end

if isfield(opt, 'epsilon')
    epsilon = opt.epsilon;
else
    epsilon = 10^-5;
end

if isfield(opt, 'max_itr')
    max_itr = opt.max_itr;
else
    max_itr = 100;
end

alpha = 1;
Z_old = W_init;
Wp = W_init;
L = L0;
fval_old = rand(1,1);
fval = rand(1,1);
itr_counter = 0;

fval_vec = [];

while abs((fval_old-fval)/fval_old)>epsilon
    itr_counter = itr_counter+1;
    fval_old = fval;
    W_old = Wp;
    [Wp,P,sval] = ComputeQP(Xtrain,Ytrain,Z_old,L,lambda);
    f = 0.5*norm(Ytrain-Xtrain*Wp,'fro')^2;
    fval = f+lambda*sval;
    Q = P+lambda*sval;
    
    while fval>Q
        fprintf('Searching step size (fval = %f, Q = %f)...\n',fval,Q);
        L = L*gamma;
        [Wp,P,sval] = ComputeQP(Xtrain,Ytrain,Z_old,L,lambda);
        f = 0.5*norm(Ytrain-Xtrain*Wp,'fro')^2;
        fval = f+lambda*sval;
        Q = P+lambda*sval;
    end

    fval_vec = [fval_vec,fval];
    
    alpha_old = alpha;
    alpha = (1+sqrt(1+4*alpha^2))/2;
    Z_old = Wp+((alpha_old-1)/alpha)*(Wp-W_old);
    
    if itr_counter>max_itr
        break;
    end
    
	fprintf('Iteration = %8d,  objective = %f\n',itr_counter, fval);
end
W = Wp;
return;

function [Wp,P,sval] = ComputeQP(X,Y,W,L,lambda)

[W1,delta_W] = ComputeGradStep(X,Y,W,L);
[U,D,V] = svd(W1,0);
D = diag(D);
D = D-(lambda/L);
idx = find(D>0);
sval = sum(D(idx));
Wp = U(:,idx)*diag(D(idx))*V(:,idx)';
P = 0.5*norm(Y-X*W,'fro')^2+trace(delta_W'*(Wp-W))+0.5*L*norm(Wp-W,'fro')^2;
return;

function [W1,delta_W] = ComputeGradStep(X,Y,W,L)

delta_W = ComputeDerivative(X,Y,W);
W1 = W-(1/L)*delta_W;
return;

function dev = ComputeDerivative(X,Y,W)

dev = X'*(X*W-Y);
return;




