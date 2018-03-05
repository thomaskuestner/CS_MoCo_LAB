function [w,fun,time,iter] = gistLogistic(X,y,lambda,theta,varargin)

% Generalized Iterative Shrinkage and Thresholding (GIST) with Logsitic Regression loss
%
% Non-convex optimization problem:
%
% min_w L(w) + \sum_i r_i(w)
%
% ================================ loss function ==========================
%
% L(w) = 1/n \sum_j log(1 + exp(-y_j*x_j'*w)) (n: number of samples)
%
% ================================ regularizer ============================
%
%  regtype = 1: Capped L1 regularizer (CapL1) (default) 
%            r_i(w) = lambda*\min(|w_i|,theta), (theta > 0, lambda >= 0)
% 
%  regtype = 2: Log Sum Penalty (LSP)
%            r_i(w) = lambda*\sum_i log(1 + |w_i|/theta), (theta > 0, lambda >= 0)
% 
%  regtype = 3: Smoothly Clipped Absolute Deviation (SCAD)
%            r_i(w) = lambda*|w_i|, if |w_i|<=lambda
%            r_i(w) = (-w_i^2 + 2*theta*lambda*|w_i| - lambda^2)/(2(theta - 1)), if lambda<=|w_i|<=theta*lambda
%            r_i(w) = 0.5*(theta + 1)*lambda^2, if |w_i| > theta*lambda, (theta > 2, lambda >= 0)
%
%  regtype = 4: Minimax Concave Penalty (MCP)
%            r_i(w) = lambda*|w_i| - 0.5*w_i^2/theta, if |w_i|<=theta*lambda
%            r_i(w) = 0.5*theta*lambda^2, if |w_i| > theta*lambda, (theta >
%            0, lambda >= 0)
%
% ============================ Input ======================================
%
% X: data matrix with each row as a sample
%
% y: label vector (+1 or -1)
%
% lambda: regularization parameter
%
% theta: theresholding parameter
%
% ======================= varargin: optional settings  ====================
%
% 'regtype': nonconvex regularization type 
%          1: CapL1 (default) 
%          2: LSP  
%          3: SCAD
%          4: MCP 
%
% 'stopcriterion': stopping criterion 
%                1: relative difference of objective functions 
%                   is less than tol (default)
%                0: relative difference of iterative weights is less
%                   than tol
%
% 'startingpoint': starting point (default: zero vector)
%
% 'tolerance': stopping tolerance (default: 1e-5)
%
% 'maxiteration': number of maximum iteration (default: 1000)
%
% 'tinitialization': initialization of t (default: 1)
%
% 'tmin': tmin parameter (default: 1e-20)
%
% 'tmax': tmax parameter (default: 1e-20)
%
% 'eta': eta factor (default: 2)
%
% 'sigma': parameter in the line search (default: 1e-5)
%
% 'nonmonotone': nonmonotone steps in the line search (default: 5)
% 
% 'stopnum': number of satisfying stopping criterion (default: 3)
%
% 'maxinneriter': number of maximum inner iteration (line search) (default: 20)
%
% ============================= Output ====================================
%
% w: output weight vector
%
% fun: a vector including all function values at each iteration
%
% time: a vector including all CPU times at each iteration
%
% iter: the number of iterative steps 
%
% =========================================================================
%
% Copyright (C) 2012-2013 Pinghua Gong
%
% For any problem, please contact Pinghua Gong via pinghuag@gmail.com
%
% Last modified on March 14, 2013.
%
% Related papers:
%
% [1] Pinghua Gong, Changshui Zhang, Zhaosong Lu, Jianhua Huang, Jieping Ye,
%     A General Iterative Shrinkage and Thresholding Algorithm for Non-convex
%     Regularized Optimization Problems. ICML 2013.
%
% ========================================================================

if nargin < 4
    error('Too few input parameters!');
end

if theta <= 0 || lambda < 0
    error('\theta must be positive and \lambda must be nonneagtive!');
end

% Parse the optional inputs.
if (mod(length(varargin), 2) ~= 0 ),
    error(['Optional Parameters passed to the function ''' mfilename ''' must be passed in pairs!']);
end

% default parameter settings
regtype = 1;
[n,d] = size(X); 
w0 = zeros(d,1);
stopcriterion = 1;
tol = 1e-5; 
maxiter = 1000;
M = 5;
sigma = 1e-5;

t = 1;
tmin = 1e-20;
tmax = 1e20;

% % t, tmin and tmax are adaptively estimated
% Linfnorm = full(max(sum(abs(X),2)))/n; L1norm = full(max(sum(abs(X),1)))/n; Lmaxnorm = full(max(max(abs(X))))/n;
% t = max([Linfnorm*Linfnorm/d,L1norm*L1norm/n,Lmaxnorm]);
% tmin = min([Linfnorm*Linfnorm/d,L1norm*L1norm/n,Lmaxnorm]);
% tmax = min([Linfnorm*L1norm,n*d*Lmaxnorm,n*Linfnorm*Linfnorm,d*L1norm*L1norm])/(1-sigma);

eta = 2;
stopnum = 3;
maxinneriter = 20;

% Optional parameter settings
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch lower(parameterName)
        case 'regtype'
            regtype = parameterValue;
            if regtype == 3 && theta <= 2 
                error('\theta must be greater than 2!');
            end
        case 'startingpoint'
            w0 = parameterValue;
        case 'stopcriterion'
            stopcriterion = parameterValue;
        case 'tolerance'
            tol = parameterValue;
        case 'maxiteration'
            maxiter = parameterValue;
        case 'nonmonotone'
            M = parameterValue;
        case 'tinitialization'
            t = parameterValue;
        case 'tmin'
            tmin = parameterValue;
        case 'tmax'
            tmax =  parameterValue;
        case 'sigma'
            sigma =  parameterValue;
        case 'eta'
            eta = parameterValue;
        case 'stopnum'
            stopnum = parameterValue;
        case 'maxinneriter'
            maxinneriter = parameterValue;
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''!']);
    end
end

w = w0; 
fun = zeros(maxiter+1,1); time = fun;
logist = zeros(n,1);

% Initial function value
Z = sparse(1:n,1:n,y,n,n)*X;
Zw = -Z*w; posind = (Zw > 0);
logist(posind) = 1 + exp(-Zw(posind));
logist(~posind) = 1 + exp(Zw(~posind));

temp = logist;
temp(posind) = 1./logist(posind);
temp(~posind) = (logist(~posind)-1)./logist(~posind);
grad =  -Z'*temp/n;

fun(1) = (sum(log(logist(~posind))) + sum(Zw(posind) + log(logist(posind))))/n + funRegC(w,d,lambda,theta,regtype);
time(1) = 0;

count = 0;
for iter = 1:maxiter
    tic;
    
    w_old = w;
    grad_old = grad;
    t = min(max(t,tmin),tmax);
    
    % line search 
    for inneriter = 1:maxinneriter
        w = proximalRegC(w_old - grad_old/t, d, lambda/t, theta,regtype);
        dw = w - w_old;
        Zw = -Z*w; posind = (Zw > 0);
        logist(posind) = 1 + exp(-Zw(posind));
        logist(~posind) = 1 + exp(Zw(~posind));
        fun(iter+1) = (sum(log(logist(~posind))) + sum(Zw(posind) + log(logist(posind))))/n + funRegC(w,d,lambda,theta,regtype);  
        if fun(iter+1) <= max(fun(max(iter-M+1,1): iter)) - 0.5*sigma*t*norm(dw)^2
            break;
        else
            t = t*eta;
        end
    end
    time(iter+1) = time(iter) + toc; 
    
    % stopping condition
    if stopcriterion
        relativediff = abs(fun(iter) - fun(iter+1))/fun(iter+1);
    else
        relativediff = norm(w - w_old)/norm(w);
    end
    if relativediff < tol
        count = count + 1;
    else
        count = 0;
    end
    if count >= stopnum
        break;
    end
    
    temp = logist;
    temp(posind) = 1./logist(posind);
    temp(~posind) = (logist(~posind)-1)./logist(~posind);
    grad =  -Z'*temp/n;
    
    % BB rule
    st = w - w_old;
    rt = grad - grad_old;
    if norm(st)/d < 1e-12 || norm(rt)/d < 1e-12
        break;
    end
    t = st'*rt/(st'*st);
    
end
fun = fun(1: min(maxiter,iter)+1);
time = time(1: min(maxiter,iter)+1);

