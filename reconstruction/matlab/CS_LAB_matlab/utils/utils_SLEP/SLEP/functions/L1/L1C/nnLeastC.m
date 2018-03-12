function [x, funVal]=nnLeastC(A, y, z, opts)
%
%%
% Function LeastC
%      Least Squares Loss with the L1 ball constraint (also known as Lasso)
%                              and the non-negative constraint
%
%% Problem
%
%  min  1/2 || A x - y||^2 + 1/2 rsL2 * ||x||_2^2
%  s.t. ||x||_1 <= z, x>=0
%
%  By default, rsL2=0.
%  When rsL2 is nonzero, this correspons the well-know elastic net.
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size mx1)
%  z -        Radius of the L1 ball (z >0)
%  opts-      Optional inputs (default value: opts=[])
%
%% Output parameters:
%  x-         Solution
%  funVal-    Function value during iterations
%
%% Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
%
% You are suggested to first read the Sparse Learning Library Manual.
%
% For any problem, please contact with Jun Liu via j.liu@asu.edu
%
% Last modified on February 18, 2010.
%
%% Related papers
%
% [1]  Jun Liu and Jieping Ye, Efficient Euclidean Projections
%      in Linear Time, ICML 2009.
%
% [2]  Jun Liu and Jieping Ye, Sparse Learning with Efficient Euclidean
%      Projections onto the L1 Ball, Technical Report ASU, 2008.
%
%% Related functions
%
%  sll_opts, initFactor, pathSolutionLeast
%  LeastR, nnLeastR, LeastC,
%  eplb
%
%%

%% Verify and initialize the parameters
%%

% Verify the number of input parameters
if (nargin <3)
    error('\n Inputs: A, y and z should be specified!\n');
elseif (nargin==3)
    opts=[];
end

% Get the size of the matrix A
[m,n]=size(A);

% Verify the length of y
if (length(y) ~=m)
    error('\n Check the length of y!\n');
end

% Verify the value of z
if (z<=0)
    error('\n z should be positive!\n');
end

% run sll_opts to set default values (flags)
opts=sll_opts(opts); 

%% Detailed initialization
%% Normalization

% Please refer to sll_opts for the definitions of mu, nu and nFlag
%
% If .nFlag =1, the input matrix A is normalized to
%                     A= ( A- repmat(mu, m,1) ) * diag(nu)^{-1}
%
% If .nFlag =2, the input matrix A is normalized to
%                     A= diag(nu)^{-1} * ( A- repmat(mu, m,1) )
%
% Such normalization is done implicitly
%     This implicit normalization is suggested for the sparse matrix
%                                    but not for the dense matrix
%

if (opts.nFlag~=0)
    if (isfield(opts,'mu'))
        mu=opts.mu;
        if(size(mu,2)~=n)
            error('\n Check the input .mu');
        end
    else
        mu=mean(A,1);
    end

    if (opts.nFlag==1)
        if (isfield(opts,'nu'))
            nu=opts.nu;
            if(size(nu,1)~=n)
                error('\n Check the input .nu!');
            end
        else
            nu=(sum(A.^2,1)/m).^(0.5); nu=nu';
        end
    else % .nFlag=2
        if (isfield(opts,'nu'))
            nu=opts.nu;
            if(size(nu,1)~=m)
                error('\n Check the input .nu!');
            end
        else
            nu=(sum(A.^2,2)/n).^(0.5);
        end
    end

    ind_zero=find(abs(nu)<= 1e-10);    nu(ind_zero)=1;
    % If some values in nu is typically small, it might be that,
    % the entries in a given row or column in A are all close to zero.
    % For numerical stability, we set the corresponding value to 1.
end

if (~issparse(A)) && (opts.nFlag~=0)
    fprintf('\n -----------------------------------------------------');
    fprintf('\n The data is not sparse or not stored in sparse format');
    fprintf('\n The code still works.');
    fprintf('\n But we suggest you to normalize the data directly,');
    fprintf('\n for achieving better efficiency.');
    fprintf('\n -----------------------------------------------------');
end

%% Starting point initialization

% compute AT y
if (opts.nFlag==0)
    ATy=A'*y;
elseif (opts.nFlag==1)
    ATy=A'*y - sum(y) * mu';  ATy=ATy./nu;
else
    invNu=y./nu;              ATy=A'*invNu-sum(invNu)*mu';
end

% L2 norm regularization
if isfield(opts,'rsL2')
    rsL2=opts.rsL2;
    if (rsL2<0)
        error('\n opts.rsL2 should be nonnegative!');
    end
else
    rsL2=0;
end

% initialize a starting point
if opts.init==2
    x=zeros(n,1);
else
    if isfield(opts,'x0')
        x=opts.x0;
        if (length(x)~=n)
            error('\n Check the input .x0');
        end
    else
        x=abs(ATy);  % if .x0 is not specified, we use ratio*ATy,
                     % where ratio is a positive value
    end
end

% compute A x
if (opts.nFlag==0)
    Ax=A* x;
elseif (opts.nFlag==1)
    invNu=x./nu; mu_invNu=mu * invNu;
    Ax=A*invNu -repmat(mu_invNu, m, 1);
else
    Ax=A*x-repmat(mu*x, m, 1);     Ax=Ax./nu;
end

if (opts.init==0) % If .init=0, we set x=ratio*x by "initFactor"
                  % Please refer to the function initFactor for detail

    x_norm=sum(abs(x)); % L1 norm of x
    x_2norm=x'*x;       % squared two norm of x
    if x_norm>=1e-6
        ratio=initFactor(x_norm, Ax, y, z,'nnLeastC', rsL2, x_2norm);
        x=ratio*x;    Ax=ratio*Ax;
    end
end

%% The main program

bFlag=0; % this flag tests whether the gradient step only changes a little

gamma=1 + rsL2;
% We assume that the maximum eigenvalue of A'A is over 1

% assign xp with x, and Axp with Ax.
% xxp=x - xp
xp=x; Axp=Ax; xxp=zeros(n,1);

%% The Armijo Goldstein line search schemes

%% if (opts.lFlag==0)

 % tp and t are used for computing the weight in forming search point
tp=0; t=1;

% lambda0 is a guess of the root in the Euclidean projection
lambda0=0;

for iterStep=1:opts.maxIter
    % --------------------------- step 1 ---------------------------
    % compute search point s based on xp and x (with alpha)
    alpha=(tp-1)/t;    s=x + alpha* xxp;

    % --------------------------- step 2 ---------------------------
    % line search for gamma and compute the new approximate solution x

    % compute the gradient (g) at s
    As=Ax + alpha* (Ax-Axp);

    % compute AT As
    if (opts.nFlag==0)
        ATAs=A'*As;
    elseif (opts.nFlag==1)
        ATAs=A'*As - sum(As) * mu';  ATAs=ATAs./nu;
    else
        invNu=As./nu;                ATAs=A'*invNu-sum(invNu)*mu';
    end

    % obtain the gradient g
    g=ATAs-ATy + rsL2 * s;

    % copy x and Ax to xp and Axp
    xp=x;    Axp=Ax;

    while (1)
        % let s walk in a step in the antigradient of s to get v
        % and project v onto the L1 ball
        v=s-g/gamma;
        p_ind=v>=0;        vv=v(p_ind); % get the positive ones
                                        % the negative ones are directly
                                        % thresholded to zero
        
        % projection (considering the constraint x>=0
        [xx, lambda, zf_step]=eplb(vv, length(vv), z, lambda0);
        lambda0=lambda;
        
        x=zeros(n,1);
        x(p_ind,1)=xx;

        v=x-s;  % the difference between the new approximate solution x
                % and the search point s

        % compute A x
        if (opts.nFlag==0)
            Ax=A* x;
        elseif (opts.nFlag==1)
            invNu=x./nu; mu_invNu=mu * invNu;
            Ax=A*invNu -repmat(mu_invNu, m, 1);
        else
            Ax=A*x-repmat(mu*x, m, 1);     Ax=Ax./nu;
        end

        Av=Ax -As;
        r_sum=v'*v; l_sum=Av'*Av;
        
        if (r_sum <=1e-20)
            bFlag=1; % this shows that, the gradient step makes little improvement
            break;
        end

        % the condition is ||Av||_2^2 <= (gamma - rsL2) * ||v||_2^2
        if(l_sum <= r_sum * (gamma-rsL2))
            break;
        else
            gamma=max(2*gamma, l_sum/r_sum + rsL2);
            % fprintf('\n gamma=%5.6f',gamma);
        end
    end

    % --------------------------- step 3 ---------------------------
    % update t and tp, and check whether converge
    tp=t; t= (1+ sqrt(4*t*t +1))/2;

    xxp=x-xp;   Axy=Ax-y;
    funVal(iterStep)=Axy' * Axy/2 + rsL2/2 * x'*x;
    
    if (bFlag)
        % fprintf('\n The program terminates as the gradient step changes the solution very small.');
        break;
    end

    switch(opts.tFlag)
        case 0
            if iterStep>=2
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <= opts.tol)
                    break;
                end
            end
        case 1
            if iterStep>=2
                if (abs( funVal(iterStep) - funVal(iterStep-1) ) <=...
                        opts.tol* funVal(iterStep-1))
                    break;
                end
            end
        case 2
            if ( funVal(iterStep)<= opts.tol)
                break;
            end
        case 3
            norm_xxp=sqrt(xxp'*xxp);
            if ( norm_xxp <=opts.tol)
                break;
            end
        case 4
            norm_xp=sqrt(xp'*xp);    norm_xxp=sqrt(xxp'*xxp);
            if ( norm_xxp <=opts.tol * max(norm_xp,1))
                break;
            end
        case 5
            if iterStep>=opts.maxIter
                break;
            end
    end
end

%% The case for .lFlag=1 is still under construction!
% We shall give the code for .lFlag=1 later.

