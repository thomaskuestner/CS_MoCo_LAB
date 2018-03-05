function [x, funVal, ValueL]=mcLeastR(A, y, z, opts)
%
%%
% Function mcLeastR:
%      Least Squares Loss for Multi-class (task) Learning
%             via the (group) L1/Lq-norm Regularization
%
%% Problem
%
%  min  1/2 || A x - y||^2 + z * sum_j ||x^j||_q
%
%  x^j denotes the j-th row of x
%
%  When y are of binary values (1 and -1), this problem is the
%  multi-class classificaiton problem, which learns the projection
%  matrix x, by sharing information across k classification tasks (with the
%  L1/Lq-norm Regularization)
%
%  The elementes in y can be any real value. In this case, this
%  problem can be regarded as the multi-task learning problem, sharing the
%  same data matrix.

%  For the case that the multi tasks have separate data
%  matrix, please refer to the functions:
%            mtLeastR and mtLogisticR.
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size mxk)
%  z -        L1/Lq norm regularization parameter (z >=0)
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
% Last modified on February 19, 2010.
%
%% Related papers
%
% [1]  Jun Liu, Shuiwang Ji, and Jieping Ye, Multi-Task Feature Learning
%      Via Efficient L2,1-Norm Minimization, UAI, 2009
%
% [2]  Jun Liu, Lei Yuan, Songcan Chen and Jieping Ye, Multi-Task Feature Learning
%      Via Efficient L2,1-Norm Minimization, Technical Report ASU, 2009.
%
%% Related functions:
%
%  sll_opts, initFactor, pathSolutionLeast
%  mcLogisticR, eppVectorR
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

if (z<0)
    error('\n z should be positive!\n');
end

opts=sll_opts(opts); % run sll_opts to set default values (flags)


% restart the program for better efficiency
%  this is a newly added function
if (~isfield(opts,'rStartNum'))
    opts.rStartNum=opts.maxIter;
else
    if (opts.rStartNum<=0)
        opts.rStartNum=opts.maxIter;
    end
end

%% Detailed initialization
%%

% Initialize q
if (~isfield(opts,'q'))
    q=2; opts.q=2;
else
    q=opts.q;
    if (q<1)
        error('\n q should be larger than 1');
    end
end

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
    ATy =A'* y;
elseif (opts.nFlag==1)
    ATy= A'* y -  mu' * sum(y, 1);    ATy=ATy./repmat(nu, 1, k);
else
    invNu=y./repmat(nu, 1, k);        ATy=A'*invNu-mu' * sum(invNu, 1);
end

% process the regularization parameter
if (opts.rFlag==0)
    lambda=z;
else % z here is the scaling factor lying in [0,1]
    if (z<0 || z>1)
        error('\n opts.rFlag=1, and z should be in [0,1]');
    end

    if q==1
        q_bar=Inf;
    elseif q>=1e6
        q_bar=1;
    else
        q_bar=q/(q-1);
    end

    lambda_max=0;
    for i=1:n
        lambda_max=max(lambda_max,...
            norm(  ATy(i,:), q_bar) );
    end
    lambda=z*lambda_max;
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
        x_norm=x_norm+ norm( x(i,:), q );
    end

    if x_norm>=1e-6
        ratio=initFactor(x_norm, Ax, y, lambda,'mcLeastR');
        x=ratio*x;    Ax=ratio*Ax;
    end
end

%% The main program

bFlag=0; % this flag tests whether the gradient step only changes a little


% The Armijo Goldstein line search scheme + accelearted gradient descent
if (opts.mFlag==0 && opts.lFlag==0)
    

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
            % and then do the L1/Lq-norm regularized projection
            v=s-g/L;

            % L1/Lq-norm regularized projection
            x=eppMatrix(v, n, k, lambda/ L, q);

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

        for i=1:n
            funVal(iterStep)=funVal(iterStep)+ lambda* norm(  x(i,:), q);
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
end

%% Reformulated problem + adpative line search

% .mFlag=1, .lFlag=1, .q=2
%  refomulate the problem as the constrained convex optimization
%  problem, and then apply adaptive line search scheme

% Original Problem:
%   min  1/2 || A x - y||^2 + z * sum_j ||x^j||_2
%
% Reformulated Problem:
%    min  1/2 || A x - y||^2 + z * t' * 1
%    s.t.   ||x^j||_2 <= t_j

% we set gamma_0 to the L that is appropriate for the starting point
% opts.x0

if (opts.mFlag==1 && opts.lFlag==1 && opts.q==2)

    L=1;
    % We assume that the maximum eigenvalue of A'A is over 1
        
    gamma=1;
    % we shall set the value of gamma = L,
    % and L is appropriate for the starting point

    xp=x; Axp=Ax;
    % store x and Ax
    xxp=zeros(n,k);
    % the difference of x and xp    
    t=sqrt(sum(x.^2, 2)); tp=t;
    % t is the upper bound of absolute value of x

    % compute AT Ax
    if (opts.nFlag==0)
        ATAx=A'*Ax;
    elseif (opts.nFlag==1)
        ATAx=A'*Ax - sum(Ax) * mu';  ATAx=ATAx./nu;
    else
        invNu=Ax./nu;                ATAx=A'*invNu-sum(invNu)*mu';
    end

    % We begin the adaptive line search in the following
    %
    % Note that, in the line search, L and beta are changing

    for iterStep=1:opts.maxIter

        ATAxp=ATAx;
        % store ATAx to ATAxp

        if (iterStep~=1)
            % compute AT Ax
            if (opts.nFlag==0)
                ATAx=A'*Ax;
            elseif (opts.nFlag==1)
                ATAx=A'*Ax - sum(Ax) * mu';  ATAx=ATAx./nu;
            else
                invNu=Ax./nu;                ATAx=A'*invNu-sum(invNu)*mu';
            end            
        end

        %--------- Line Search for L begins        
        while (1) 
            if (iterStep~=1)
                alpha= (-gamma+ sqrt(gamma*gamma + 4* L * gamma)) / (2*L);
                beta= (gamma - gamma* alphap) / (alphap * gamma + alphap* L * alpha);
                % beta is the coefficient for generating search point s

                s=x + beta* xxp;   s_t= t + beta * (t -tp);
                As=Ax + beta* (Ax-Axp);
                ATAs=ATAx + beta * (ATAx- ATAxp);
                % compute the search point s, A * s, and A' * A * s
            else
                alpha= (-1+ sqrt(5)) / 2;
                beta=0; 
                
                s=x; s_t=t; As=Ax; ATAs=ATAx;
            end

            g=ATAs-ATy;
            % compute the gradient g

            % let s walk in a step in the antigradient of s
            u=s-g / L;
            v= s_t - lambda / L;

            % projection
            [xnew, tnew]=ep21R(u,v,n,k);

            v=xnew-s;  % the difference between the new approximate solution x
            % and the search point s
            v_t=tnew-s_t;

            % compute A * xnew
            if (opts.nFlag==0)
                Axnew=A* xnew;
            elseif (opts.nFlag==1)
                invNu=xnew./nu; mu_invNu=mu * invNu;
                Axnew=A*invNu -repmat(mu_invNu, m, 1);
            else
                Axnew=A*xnew-repmat(mu*xnew, m, 1);     Axnew=Axnew./nu;
            end

            Av=Axnew -As;
            r_sum=norm(v,'fro')^2 + v_t'*v_t; l_sum=norm(Av,'fro')^2;
            
            if (r_sum <=1e-20)
                bFlag=1; % this shows that, the gradient step makes little improvement
                break;
            end            

            % the condition is ||Av||_2^2
            %                       <= L * (||v||_2^2 + ||v_t|| _2^2 )
            if(l_sum <= r_sum * L)
                break;
            else
                L=max(2*L, l_sum/r_sum);
                % fprintf('\n L=%5.6f',L);
            end
        end
        %--------- Line Search for L ends
            
        gamma=L* alpha* alpha;    alphap=alpha;
        % update gamma, and alphap

        ValueL(iterStep)=L;
        % store values for L

        tao=L * r_sum / l_sum;
        if (tao >=5)
            L=L*0.8;
        end
        % decrease the value of L

        xp=x;   x=xnew;  xxp=x-xp; 
        Axp=Ax; Ax=Axnew;
        % update x and Ax with xnew and Axnew
        tp=t; t=tnew;
        % update tp and t
        
        Axy=Ax-y;        
        funVal(iterStep)=norm(Axy,'fro')^2/2 + lambda * sum(t);
        % compute function value

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
    end
end

%% Reformulated problem + Nemirovski's line search

% .mFlag=1, .lFlag=0, .q=2
%  refomulate the problem as the constrained convex optimization
%  problem, and then apply Armijo Goldstein line search scheme

% Original Problem:
%   min  1/2 || A x - y||^2 + z * sum_j ||x^j||_2
%
% Reformulated Problem:
%    min  1/2 || A x - y||^2 + z * t' * 1
%    s.t.   ||x^j||_2 <= t_j

if (opts.mFlag==1 && opts.lFlag==0 && opts.q==2)

    L=1;
    % We assume that the maximum eigenvalue of A'A is over 1

    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,k);

    alphap=0; alpha=1;
    
    t=sqrt(sum(x.^2, 2)); tp=t;
    % t is the upper bound of absolute value of x

    for iterStep=1:opts.maxIter
        % --------------------------- step 1 ---------------------------
        % compute search point s based on xp and x (with beta)
        beta=(alphap-1)/alpha;    s=x + beta* xxp; s_t= t + beta * (t-tp);

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
        xp=x;    Axp=Ax; tp=t;

        while (1)
            % let (s, s_t) walk in a step in the antigradient of s to get
            % (u,v), and then do the L1/Lq-norm regularized projection
            u=s-g/L;     v= s_t - lambda / L;

            % projection
            [x, t]=ep21R(u,v,n,k);

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
            r_sum=norm(v,'fro')^2 + norm(t-s_t)^2; l_sum=norm(Av,'fro')^2;
            
            if (r_sum <=1e-20)
                bFlag=1; % this shows that, the gradient step makes little improvement
                break;
            end            

            % the condition is ||Av||_2^2
            %                       <= L * (||v||_2^2 + ||t-s_t|| _2^2 )
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
        funVal(iterStep)=norm(Axy,'fro')^2/2 +  lambda * sum(t);

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
    end
end

%%
if(opts.mFlag==0 && opts.lFlag==1)
    error('\n The function does not support opts.mFlag=0 & opts.lFlag=1!');
end