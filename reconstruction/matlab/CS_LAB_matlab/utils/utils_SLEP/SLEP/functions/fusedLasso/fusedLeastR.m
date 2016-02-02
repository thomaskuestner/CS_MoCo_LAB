function [x, funVal, ValueL]=fusedLeastR(A, y, lambda, opts)
%
%%
% Function fusedLeastR
%      Least Squares Loss with the Fused Lasso Penalty
%
%% Problem
%
%  min  1/2 || A x - y||^2  + lambda * ||x||_1 +   lambda_2 ||Rx||_1
%                               %% + 1/2 rsL2 * ||x||_2^2
%
%  %%By default, rsL2=0.
%  %%When rsL2 is nonzero, this correspons the well-know elastic net.
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size mx1)
%  lambda -        L_1 norm regularization parameter (lambda >=0)
%  opts-      Optional inputs (default value: opts=[])
%
%% Output parameters:
%  x-         Solution
%  funVal-    Function value during iterations
%
%% Related papers
%
% [1]  Jun Liu, Lei Yuan, and Jieping Ye, An Efficient Algorithm for 
%      a Class of Fused Lasso Problems, KDD, 2010.
%
%% Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
%
% You are suggested to first read the Manual.
%
% For any problem, please contact with Jun Liu via j.liu@asu.edu
%
% Last modified 14 November 2009.
%
% Related functions:
%  sll_opts, initFactor, LeastC
%%

%% Verify and initialize the parameters
%%

% Verify the number of input parameters
if (nargin <3)
    error('\n Inputs: A, y and lambda should be specified!\n');
elseif (nargin==3)
    opts=[];
end

% Get the size of the matrix A
[m,n]=size(A);

% Verify the length of y
if (length(y) ~=m)
    error('\n Check the length of y!\n');
end

% Verify the value of lambda
if (lambda<0)
    error('\n lambda should be nonnegative!\n');
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

% process the regularization parameter

% % % L2 norm regularization
% % if isfield(opts,'rsL2')
% %     rsL2=opts.rsL2;
% %     if (rsL2<0)
% %         error('\n opts.rsL2 should be nonnegative!');
% %     end
% % else
% %     rsL2=0;
% % end

% L1 norm regularization
if (opts.rFlag==0)
    if ~isfield(opts,'fusedPenalty')
        lambda2=0;
    else
        lambda2=opts.fusedPenalty;
    end
else % lambda here is the scaling factor lying in [0,1]
    if (lambda<0 || lambda>1)
        error('\n opts.rFlag=1, and lambda should be in [0,1]');
    end

    lambda_max=max(abs(ATy));
    lambda=lambda*lambda_max;

    if ~isfield(opts,'fusedPenalty')
        lambda2=0;
    else
        lambda2=lambda_max * opts.fusedPenalty;
    end

    %%rsL2=rsL2*lambda_max; % the input rsL2 is a ratio of lambda_max
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
        x=ATy;  % if .x0 is not specified, we use ratio*ATy,
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

    x_norm=sum(abs(x)); x_2norm=x'*x;
    if x_norm>=1e-6
        %ratio=initFactor(x_norm, Ax, y, lambda,'LeastR', rsL2, x_2norm);
        ratio=initFactor(x_norm, Ax, y, lambda,'LeastR', 0, x_2norm);
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

z0=zeros(n-1,1);
% z0 is the starting point for flsa

bFlag=0; % this flag tests whether the gradient step only changes a little

L=1;
%%L=1 + rsL2;
% We assume that the maximum eigenvalue of A'A is over 1

%% The Armijo Goldstein line search scheme + accelearted gradient descent
if (opts.lFlag==0)

    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,1);

    % alphap and alpha are used for computing the weight in forming search point
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
            ATAs=A'*As - sum(As) * mu';  ATAs=ATAs./nu;
        else
            invNu=As./nu;                ATAs=A'*invNu-sum(invNu)*mu';
        end

        % obtain the gradient g
        %%g=ATAs-ATy + rsL2 * s;
        g=ATAs-ATy;

        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;

        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and then do the l1-norm regularized projection
            v=s-g/L;

            [x, z, infor]=flsa(v, z0,...
                lambda / L, lambda2 / L, n,...
                1000, 1e-8, 1, 6);
            z0=z;

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

            % the condition is ||Av||_2^2 <= (L - rsL2) * ||v||_2^2
            %%if(l_sum <= r_sum * (L-rsL2))
            if(l_sum <= r_sum * L)
                break;
            else
                %%L=max(2*L, l_sum/r_sum + rsL2);
                L=max(2*L, l_sum/r_sum);
                % fprintf('\n L=%5.6f',L);
            end
        end

        %ValueL(iterStep)=L;
        ValueL(iterStep)=infor(2);

        % --------------------------- step 3 ---------------------------
        % update alpha and alphap, and check whether converge
        alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;

        xxp=x-xp;   Axy=Ax-y;
        %%funVal(iterStep)=Axy'* Axy/2 + rsL2/2 * x'*x + sum(abs(x)) * lambda + lambda2 * sum( abs( x(2:n) -x(1:(n-1)) ));
        funVal(iterStep)=Axy'* Axy/2  + sum(abs(x)) * lambda + lambda2 * sum( abs( x(2:n) -x(1:(n-1)) ));
        
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
        
        % restart the program every opts.rStartNum
        if (~mod(iterStep, opts.rStartNum))
            alphap=0; alpha=1;
            xp=x; Axp=Ax; xxp=zeros(n,1); L =1;
        end
    end
end


%% the line search scheme by Nesterov
if(opts.lFlag==1)
    
    % compute AT Ax
    if (opts.nFlag==0)
        ATAx=A'*Ax;
    elseif (opts.nFlag==1)
        ATAx=A'*Ax - sum(Ax) * mu';  ATAx=ATAx./nu;
    else
        invNu=As./nu;                ATAx=A'*invNu-sum(invNu)*mu';
    end
    g_x=ATAx-ATy;  % g_x is the gradient at x
        
    % initialization
    v=x;              % v is the minimizer of \psi_i(x)
    g_v=g_x;       % g_v is the gradient at v
    sum_g=x;      % sum_g contains the summation of x0 - \sum_{i >=1} g_{x_i}
    sum_a=0;      % a is a nonnegative tuning parameter, and suma is the summation of all the
                       % a's during the iteration                       
    
    for iterStep=1:opts.maxIter
        while (1)
            % compute a, which is the solution to the quadratic function
            a= ( 1+sqrt(1+2*L*sum_a) ) / L;

            % compute the search point s
            s=sum_a/(sum_a+a) * x + a /(sum_a+a) * v;

            % the gradient at the search point s
            g_s=sum_a/(sum_a+a) * g_x + a /(sum_a+a) * g_v;

            % compute the new soltuion xnew
            u=s-g_s/L;      % u is a gradient based on s

            % obtain xnew by projection
            z0= z0 / L;
             [xnew, z, infor]=flsa(u, z0,...
                lambda / L, lambda2 / L, n,...
                1000, 1e-8, 1, 6);
            z0=z * L;

            % compute A xnew
            if (opts.nFlag==0)
                Axnew=A* xnew;
            elseif (opts.nFlag==1)
                invNu=xnew./nu; mu_invNu=mu * invNu;
                Axnew=A*invNu -repmat(mu_invNu, m, 1);
            else
                Axnew=A*xnew-repmat(mu*xnew, m, 1);     Axnew=Axnew./nu;
            end
            
            Ax=Axnew;
            % compute AT Ax
            if (opts.nFlag==0)
                ATAx=A'*Ax;
            elseif (opts.nFlag==1)
                ATAx=A'*Ax - sum(Ax) * mu';  ATAx=ATAx./nu;
            else
                invNu=Ax./nu;                ATAx=A'*invNu-sum(invNu)*mu';
            end
            g_xnew=ATAx-ATy;  % g_xnew is the gradient at xnew

            % test whether L is appropriate
            % cG denotes the composite gradient
            cG=L*(s- xnew) + g_xnew - g_s;
            val_left=cG' * (s- xnew); val_right=cG'*cG;
            if (  val_left * L >= val_right )

                ValueL(iterStep)=L;
                
                L=L/1.1;

                break;
            else
                L=max(2*L);
                %fprintf('\n L=%5.2f',L);
            end
        end

        % compute v as the minizer of \psi_k (x)
        sum_g=sum_g - g_xnew * a;
        sum_a=sum_a+a;

        % obtain v by projection
        z0=z0 * sum_a;
        [v, z, infor]=flsa(sum_g, z0,...
            lambda * sum_a, lambda2 * sum_a, n,...
            1000, 1e-8, 1, 6);
        z0=z / sum_a;
        
        % compute Av
        if (opts.nFlag==0)
            Av=A* v;
        elseif (opts.nFlag==1)
            invNu=v./nu; mu_invNu=mu * invNu;
            Av=A*invNu -repmat(mu_invNu, m, 1);
        else
            Av=A*v-repmat(mu*v, m, 1);     Av=Av./nu;
        end
        
        % compute AT Av
        if (opts.nFlag==0)
            g_v=A'*Av;
        elseif (opts.nFlag==1)
            g_v=A'*Av - sum(Av) * mu';  g_v=g_v./nu;
        else
            invNu=Av./nu;                g_v=A'*invNu-sum(invNu)*mu';
        end
        g_v=g_v-ATy;  % g_x is the gradient at x
            
        xxp=xnew-x;   norm_xp=norm(x,'fro');    norm_xxp=norm(xxp,'fro');
        x=xnew;        % Ax has been set to Axnew
        g_x=g_xnew;
        
        % compute the objective function value
        Axy=Ax-y;
        funVal(iterStep)=Axy'* Axy/2  + sum(abs(x)) * lambda + lambda2 * sum( abs( x(2:n) -x(1:(n-1)) ));

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
       
       if (~mod(iterStep, opts.rStartNum))
           v=x;              % v is the minimizer of \psi_i(x)
           g_v=g_x;       % g_v is the gradient at v
           sum_g=x;      % sum_g contains the summation of x0 - \sum_{i >=1} g_{x_i}
           sum_a=0;      % a is a nonnegative tuning parameter, and suma is the summation of all the
                              % a's during the iteration
       end
    end
end