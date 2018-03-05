function Theta=sparseInverseCovariance(S, lambda, opts)
%
% Function sparseInverseCovariance
%   sparse inverse covariance estimation
%
%% Problem
%
%   max log( det( Theta ) ) - < S, Theta> - lambda * ||Theta||_1
%
%% Input parameters:
%
% S-          Empirical covariance matrix
% lambda-     regularization parameter
% opts-       optional inputs
%
%% Output parameter:
% Theta-      sparse inverse covariance
%
%% Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
%
% You are suggested to first read the Manual.
%
% For any problem, please contact with Jun Liu via j.liu@asu.edu
%
% Composed on September 18, 2008.
%
%% Related papers
%
% The implementation is based on the following paper:
%
% [1]  Jerome Friedman, Trevor Hastie, and Robert Tibshirani,
%      Sparse inverse covariance estimation with the graphical lasso, 2007
%
% This function has been used in the following paper:
%
% [2] Shuai Huang, Jing Li, Liang Sun, Jun Liu, Teresa Wu,
%     Kewei Chen, Adam Fleisher, Eric Reiman, and Jieping Ye,
%     Learning Brain Connectivity of Alzheimer's Disease 
%     from Neuroimaging Data, NIPS, 2009
%
%% Related functions
%
%  invCov.c
%

% Verify the number of input parameters
if (nargin <2)
    error('\n Inputs: S and lambda should be specified!\n');
elseif (nargin==3)
    opts=[];
end

if ~isfield(opts,'LassoMaxIter')
    opts.LassoMaxIter=1000;
end
if ~isfield(opts,'fGap')
    opts.fGap=1e-6;
end
if ~isfield(opts,'xGap')
    opts.xGap=1e-4;
end

if ~isfield(opts,'maxIter')
    opts.maxIter=1000;
end
if ~isfield(opts,'xtol')
    opts.xtol=1e-4;
end


sum_S=sum(diag(S));
n=size(S,1);

% initialize Theta and a working variable W
Theta=zeros(n,n);
W=zeros(n,n);

% call the C function
invCov(Theta, W, S, lambda, sum_S, n, opts.LassoMaxIter, opts.fGap, opts.xGap, opts.maxIter, opts.xtol);

%% for the Lasso (by CD)
% LassoMaxIter: 
% fGap
% xGap

% Note:
% the Lasso terminates either
% the change of the function value is less than fGap
% or the change of the solution in terms of L1 norm is less than xGap
% or the maximal iteration maxIter has achieved


%% for the outer loop
% maxIter
% xtol

% Note
% The outer loop terminates either the difference between ajacent solution in terms of L1 norm is less than xtol,
% or the maximal iterations has achieved