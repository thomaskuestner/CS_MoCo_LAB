% revised on September 21, 2010 to make parameters consistent

% September 25, 2009
% written by Shuiwang Ji and Jieping Ye

% This function implements the accelerated gradient algorithm for matrix
% classification described in Ji and Ye (ICML 2009). The original
% formulation is described in Tomioka and Aihara (ICML 2007) in which
% each input data sample consists of a matrix and the outputs are binary
% classes. Note that the current version of code only works for binary-class
% problems.

% References:
%Ji, S. and Ye, J. 2009. An accelerated gradient method for trace norm minimization. 
%In Proceedings of the 26th Annual international Conference on Machine Learning 
%(Montreal, Quebec, Canada, June 14 - 18, 2009). ICML '09, vol. 382. ACM, New York, 
%NY, 457-464.

%Tomioka, R. and Aihara, K. 2007. Classifying matrices with a spectral regularization. 
%In Proceedings of the 24th international Conference on Machine Learning 
%(Corvalis, Oregon, June 20 - 24, 2007). Z. Ghahramani, Ed. ICML '07, vol. 
%227. ACM, New York, NY, 895-902.


%[W,b,fval_vec,itr_counter] =
%accel_grad_mc(A,Y,lambda,opt)

% required inputs:
% A: C x C x N array where each data point is a CxC matrix and N is the sample size
% Y: N x 1 vector of class labels
% lambda: regularization parameter

% optional inputs:
% opt.L0: Initial guess for the Lipschitz constant
% opt.gamma: the multiplicative factor for Lipschitz constant
% opt.W_init: initial weight matrix
% opt.b_init: initial value for bias
% opt.epsilon: precision for termination
% opt.max_itr: maximum number of iterations
% opt.loss_prim: value of the primal objective function for termination

% outputs:
% W: the computed weight matrix
% b: the computed bias
% fval_vec: a vector for the sequence of function values
% itr_counter: number of iterations executed

function [W,b,fval_vec,itr_counter] = accel_grad_mc(A,Y,lambda,opt)

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
    W_init = zeros(size(Xtrain,2),size(Xtrain,1));
end

if isfield(opt, 'b_init')
    b_init = opt.b_init;
else
    b_init = 0;
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

if isfield(opt, 'loss_prim')
    loss_prim = opt.loss_prim;
else
    loss_prim = -1;
end
    
fval_vec = [];

W_old = W_init;
L = L0;
fval_old = rand(1,1);
c_init = b_init;
fval = loss_prim+100;
itr_counter = 0;
Z_old = W_old;
alpha = 1;
loss_prim = loss_prim+10^-2;

while (fval>loss_prim)&&(abs((fval_old-fval)/fval_old)>epsilon)&&(itr_counter<max_itr)
    itr_counter = itr_counter+1;
    fval_old = fval;
    [Wp,b1,P,sval] = ComputeQP(Xtrain,Ytrain,Z_old,c_init,L,lambda);
    f = ComputeFun(Xtrain,Ytrain,Wp,b1);
    fval = f+lambda*sval;
    Q = P+lambda*sval;
    
    while fval>Q
        fprintf('Searching step size (fval = %f, Q = %f)...\n',fval,Q);
 %       fprintf('norm W = %f\n',norm(Wp,'fro'));
        L = L*gamma;
        [Wp,b1,P,sval] = ComputeQP(Xtrain,Ytrain,Z_old,c_init,L,lambda);
        f = ComputeFun(Xtrain,Ytrain,Wp,b1);
        fval = f+lambda*sval;
        Q = P+lambda*sval;
    end

    fval_vec = [fval_vec,fval];
    
    alpha_old = alpha;
    alpha = (1+sqrt(1+4*alpha_old^2))/2;
    Z_old = Wp+((alpha_old-1)/alpha)*(Wp-W_old);
    c_init = b1+((alpha_old-1)/alpha)*(b1-b_init);
    b_init = b1;
    W_old = Wp;
    
    if mod(itr_counter,50)==0
        fprintf('Iteration = %8d,  objective = %f\n',itr_counter, fval);
    end
end
b = b1;
W = Wp;
return;

function [Wp,b1,P,sval] = ComputeQP(X,Y,W,b,L,lambda)

[W1,b1,delta_W,delta_b,f] = ComputeGradStep(X,Y,W,L,b);

[U,D,V] = svd(W1,0);

%[U,D] = eig(W1);

D = diag(D);
D = D-(lambda/L);
idx = find(D>0);
sval = sum(D(idx));

Wp = U(:,idx)*diag(D(idx))*V(:,idx)';

%Wp = U(:,idx)*diag(D(idx))*U(:,idx)';

%P = f+trace(delta_W'*(Wp-W))+delta_b*(b1-b)+0.5*L*(norm(Wp-W,'fro')^2+(b1-b)^2);

%disp('Compute P');
P = f+ComputeProdTrace(delta_W,(Wp-W))+delta_b*(b1-b)+0.5*L*(norm(Wp-W,'fro')^2+(b1-b)^2);

return;

function [W1,b1,delta_W,delta_b,f] = ComputeGradStep(X,Y,W,L,b)

[delta_W,delta_b,f] = ComputeDerivative(X,Y,W,b);
W1 = W-(1/L)*delta_W;
b1 = b-(1/L)*delta_b;
return;

function [dev,delta_b,f] = ComputeDerivative(X,Y,W,b)

num_samples = size(X,3);
tmp = zeros(num_samples,1);

%disp('Compute derivative');
for i = 1:num_samples
    %tmp(i) = exp(-Y(i)*(trace(W'*X(:,:,i))+b));
    tmp(i) = exp(-Y(i)*(ComputeProdTrace(W,X(:,:,i))+b));
end

dev = zeros(size(W,1),size(W,2));
delta_b = 0;
f = 0;

for i = 1:num_samples
    tt = 1+tmp(i);
    f = f+log(tt);
    dev = dev-(Y(i)*X(:,:,i)*tmp(i))/tt;
    delta_b = delta_b-Y(i)*tmp(i)/tt;
end
return;

function f = ComputeFun(X,Y,W,b)

%disp('Compute function value');
num_samples = size(X,3);
f = 0;
for i = 1:num_samples
    %f = f+log(1+exp(-Y(i)*(trace(W'*X(:,:,i))+b)));
    f = f+log(1+exp(-Y(i)*(ComputeProdTrace(W,X(:,:,i))+b)));
    
end
return;

function tval = ComputeProdTrace(A,B)

%disp('Compute trace');
%tval = 0;
%for i = 1:size(A,2)
%    tval = tval+A(:,i)'*B(:,i);
%end

a = reshape(A, [size(A,1)*size(A,2),1]);
b = reshape(B, [size(B,1)*size(B,2),1]);
tval = a'*b;
