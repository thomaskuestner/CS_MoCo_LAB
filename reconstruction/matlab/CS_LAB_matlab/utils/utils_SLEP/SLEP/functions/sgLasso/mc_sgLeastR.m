function [x, funVal, ValueL]=mc_sgLeastR(A, y, z, opts)
%
%%
% Function mc_cgLassoLeast:
%      Least Squares Loss for Multi-class (task) Learning
%             via the sparse group Lasso penalty
%
%% Problem
%
%  min  1/2 || A x - y||^2 + lambda_1 * \|x\|_1 + \lambda_2 *\|x\|_{q,1}
%
%  q=2, infty
%
%  When y are of binary values (1 and -1), this problem is the
%  multi-class classificaiton problem, which learns the projection
%  matrix x, by sharing information across k classification tasks (with the
%  composite group Lasso penalty)
%
%  The elementes in y can be any real value. In this case, this
%  problem can be regarded as the multi-task learning problem, sharing the
%  same data matrix.
%
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size mxk)
%  z -        cgLasso penalty [lambda1, lambda2]
%  opts-      Optional inputs (default value: opts=[])
%
%% Output parameters:
%  x-         Solution (of size n x k)
%  funVal-    Function value during iterations
%
%% Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
%
% You are suggested to first read the Manual.
%
% For any problem, please contact with Jun Liu via j.liu@asu.edu
%
% Last modified 1 December 2009.
%
%%
%% Related papers
%
% [1] Jun Liu and Jieping Ye, Moreau-Yosida Regularization for 
%     Grouped Tree Structure Learning, NIPS 2010
%
%% Related functions:
%
%  sll_opts
%
%%

%% Verify and initialize the parameters
%%
if (nargin <3)
    error('\n Inputs: A, y and z should be specified!\n');
end

[m,n]=size(A);

k=size(y,2);

if (size(y,1) ~=m)
    error('\n Check the length of y!\n');
end

if (length(z)~=2)
    error('\n z should contain two nonnegative parameters!\n');
end

lambda1=z(1);
lambda2=z(2);

opts=sll_opts(opts); % run sll_opts to set default values (flags)

%% Detailed initialization
%%

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

% Initialize q
if (~isfield(opts,'q'))
    opts.q=2;
else
    if (opts.q~=0 && opts.q~=2)
        error('\n opts.q should be either 0 or 2');
        % 0 stands for inf
        % 2 stands for 2
    end
end

%% Starting point initialization

% compute AT y
if (opts.nFlag==0)
    ATy =A'* y;
elseif (opts.nFlag==1)
    ATy= A'* y -  mu' * sum(y, 1);    ATy=ATy./repmat(nu, 1, k);
else
    invNu=y./repmat(nu, 1, k);        ATy=A'*invNu-mu' * sum(invNu, 1);
end

% process the regularization parameter
if (opts.rFlag~=0) % lambda1 and lambda2 are the scaling factor lying in [0,1]
    if ( lambda1<0 || lambda1>1 ||  lambda2<0 || lambda2>1)
        error('\n opts.rFlag=1, and lambda1 and lambda2 should be in [0,1]');
    end

    temp=abs(ATy);
    lambda1_max=max(max(abs(temp)));
    lambda1=lambda1*lambda1_max;
    
    temp=max(temp-lambda1,0);

    if opts.q==0
        lambda2_max=max( sum(temp,2) );
    else
        lambda2_max=sqrt( max( sum(temp.^2,2) ) );
    end
    
    lambda2=lambda2*lambda2_max;
    
    %[lambda1_max,lambda2_max, lambda1,lambda2]
end

% initialize a starting point
if opts.init==2
    x=zeros(n,k);
else
    if isfield(opts,'x0')
        x=opts.x0;
        if ( size(x,1)~=n || size(x,2)~=k )
            error('\n Check the input .x0');
        end
    else
        x=ATy;  % if .x0 is not specified, we use ratio*ATy,
        % where ratio is a positive value
    end
end

% compute Ax= A * x
if (opts.nFlag==0)
    Ax=A* x;
elseif (opts.nFlag==1)
    invNu=x./repmat(nu, 1, k);  mu_invNu=mu * invNu;
    Ax=A*invNu -repmat(mu_invNu, m, 1);
else
    Ax=A*x-repmat(mu*x, m, 1);     Ax=Ax./repmat(nu, 1, k);
end

if (opts.init==0) % If .init=0, we set x=ratio*x by "initFactor"
    % Please refer to the function initFactor for detail

    x_norm=0;
    for i=1:n
        x_norm=x_norm+ norm( x(i,:), opts.q );
    end

    if x_norm>=1e-6
        ratio=initFactor(x_norm, Ax, y, lambda,'mcLeastR');
        x=ratio*x;    Ax=ratio*Ax;
    end
end

% restart the program for better efficiency
%  this is a newly added function
if (~isfield(opts,'rStartNum'))
    opts.rStartNum=opts.maxIter;
else
    if (opts.rStartNum<=0)
        opts.rStartNum=opts.maxIter;
    end
end

%% The main program

bFlag=0; % this flag tests whether the gradient step only changes a little


% The Armijo Goldstein line search scheme + accelearted gradient descent
if (opts.lFlag==0)

    L=1;
    % We assume that the maximum eigenvalue of A'A is over 1

    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,k);

    alphap=0; alpha=1;

    for iterStep=1:opts.maxIter
        % --------------------------- step 1 ---------------------------
        % compute search point s based on xp and x (with beta)
        beta=(alphap-1)/alpha;    s=x + beta* xxp;

        % --------------------------- step 2 ---------------------------
        % line search for L and compute the new approximate solution x

        % compute the gradient (g) at s
        As=Ax + beta* (Ax-Axp);

        % compute AT As
        if (opts.nFlag==0)
            ATAs=A'*As;
        elseif (opts.nFlag==1)
            ATAs= A'* As -  mu' * sum(As, 1);  ATAs=ATAs./repmat(nu, 1, k);
        else
            invNu=As./repmat(nu, 1, k);        ATAs=A'*invNu-mu' * sum(invNu, 1);
        end

        % obtain the gradient g
        g=ATAs-ATy;

        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;

        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and then do the projection
            v=s-g/L;


            % projection
            [x,normx]=epsgLasso(v, n, k, lambda1 / L, lambda2 / L, opts.q);
            
            v=x-s;  % the difference between the new approximate solution x
            % and the search point s

            % compute Ax= A * x
            if (opts.nFlag==0)
                Ax=A* x;
            elseif (opts.nFlag==1)
                invNu=x./repmat(nu, 1, k);  mu_invNu=mu * invNu;
                Ax=A*invNu -repmat(mu_invNu, m, 1);
            else
                Ax=A*x-repmat(mu*x, m, 1);     Ax=Ax./repmat(nu, 1, k);
            end

            Av=Ax -As;
            r_sum=norm(v,'fro')^2; l_sum=norm(Av,'fro')^2;
            
            if (r_sum <=1e-20)
                bFlag=1; % this shows that, the gradient step makes little improvement
                break;
            end
            
            % the condition is ||Av||_2^2 <= L * ||v||_2^2
            if(l_sum <= r_sum * L)
                break;
            else
                L=max(2*L, l_sum/r_sum);
                % fprintf('\n L=%5.6f',L);
            end
        end
        
        ValueL(iterStep)=L;

        % --------------------------- step 3 ---------------------------
        % update alpha and alphap, and check whether converge
        alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;

        xxp=x-xp;   Axy=Ax-y;
        funVal(iterStep)=norm(Axy,'fro')^2/2;

       if (opts.q==2)
           funVal(iterStep)=funVal(iterStep)+...
               lambda1* normx(1) +...
               lambda2* normx(2);
       else
           funVal(iterStep)=funVal(iterStep)+...
               lambda1* sum(abs( x(:) ) ) +...
               lambda2* sum( max(abs(x), [], 2) );
       end
        
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
                norm_xxp=norm(xxp,'fro');
                if ( norm_xxp <=opts.tol)
                    break;
                end
            case 4
                norm_xp=norm(xp,'fro');    norm_xxp=norm(xxp,'fro');
                if ( norm_xxp <=opts.tol * max(norm_xp,1))
                    break;
                end
            case 5
                if iterStep>=opts.maxIter
                    break;
                end
        end
        
        % restart the program every opts.rStartNum
        if (~mod(iterStep, opts.rStartNum))
            alphap=0; alpha=1;
            xp=x; Axp=Ax; xxp=zeros(n,k); L =1;
        end
    end
else
    error('\n The current code is only applicable to the case of opts.lFlag=0');
end
