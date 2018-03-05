function [x, funVal, ValueL]=orderLeast(A, y, lambda, opts)
%
%%
% Function LeastC
%      Least Squares Loss with the ordered tree constraints
%
%% Problem
%
%  min  1/2 || A x - y||^2 + lambda * \|x\|
%  s.t. x satisfy the ordered tree structure (non-negative max-heap)
%
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size mx1)
%  lambda -   regularization paramter
%  opts-      Optional inputs (default value: opts=[])
%
%% Output parameters:
%  x-         Solution
%  funVal-    Function value during iterations
%
%% Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
%
% You are suggested to first read the Manual.
%
% For any problem, please contact with Jun Liu via j.liu@asu.edu
%
% Last modified on December 26, 2011.
%
%% Related papers
%
%
%% Related functions
%
%  sll_opts, initFactor, pathSolutionLeast,
%  LeastR, nnLeastR, nnLeastC,
%  eplb
%
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
if (lambda<=0)
    error('\n lambda should be positive!\n');
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

% set the regularization parameter
if isfield(opts,'rFlag')    
    if (opts.rFlag==1)
        if (lambda >1)
            error('lambda should be within [0,1]');
        end
        
        lambda_max=max(abs(ATy));
        
        lambda=lambda_max*lambda;
    end
end

% default: sequential tree
if ~isfield(opts,'treeFlag')
    opts.treeFlag=1;
else
    if (opts.treeFlag<1 || opts.treeFlag>4)
        error('opts.treeFlag is among 1,2,3,4')
    end
    
    if (opts.treeFlag==4)
        if ~isfield(opts,'FileName')
            error('\n Error! You need to specify the txt for constructing the tree');
        end
        if ~isfield(opts,'rootNum')
            error('\n Error! You need to specify the rootNum for constructing the tree');
        end
        
        FileName=opts.FileName;
        rootNum=opts.rootNum;
    end
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

    x_norm=sum(abs(x)); % L1 norm of x
    x_2norm=x'*x;       % squared two norm of x
    if x_norm>=1e-6
        ratio=initFactor(x_norm, Ax, y, lambda,'LeastR', 0, x_2norm);
        x=ratio*x;    Ax=ratio*Ax;
    end
end

%% The main program

bFlag=0; % this flag tests whether the gradient step only changes a little

%% The Armijo Goldstein line search schemes
if (opts.lFlag==0)

    L=1;
    % We assume that the maximum eigenvalue of A'A is over 1

    % assign xp with x, and Axp with Ax.
    % xxp=x - xp
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
        g=ATAs-ATy + lambda * ones(n,1);

        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;

        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and project v onto the L1 ball
            v=s-g/L;

            % projection, we need to revise this part
            switch(opts.treeFlag)
                case 1
                    x=sequence_bottomup(v,n);
                case 2
                    x=orderTreeBinary(v,n);
                case 3
                    x=orderTreeDepth1(v,n);
                case 4
                    x=orderTree(FileName, v, rootNum, n);
                otherwise
                    x=sequence_bottomup(v,n);
            end
            
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

            % the condition is ||Av||_2^2 <= (L) * ||v||_2^2
            if(l_sum <= r_sum * (L))
                break;
            else
                L=max(2*L, l_sum/r_sum );
                % fprintf('\n L=%5.6f',L);
            end
        end

        ValueL(iterStep)=L;

        % --------------------------- step 3 ---------------------------
        % update alpha and alphap, and check whether converge
        alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;

        xxp=x-xp;   Axy=Ax-y;
        funVal(iterStep)=Axy' * Axy/2 + lambda * sum(x);
        
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
    end % end of .lFlag=0
end


%% Adaptive Line Search

% we set gamma_0 to the L that is appropriate for the starting point
% opts.x0

if (opts.lFlag==1)

    L=1 ;
    % We assume that the maximum eigenvalue of A'A is over 1

    lambda0=0;
    % lambda0 is a guess of the root in the Euclidean projection

    gamma=1;
    % we shall set the value of gamma = L,
    % and L is appropriate for the starting point

    xp=x; Axp=Ax;
    % store x and Ax
    xxp=zeros(n,1);
    % the difference of x and xp

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

                s=x + beta* xxp;
                As=Ax + beta* (Ax-Axp);
                ATAs=ATAx + beta* (ATAx-ATAxp);
                % compute the search point s, A * s, and A' * A * s
            else
                alpha= (-1+ sqrt(5)) / 2;
                beta=0; s=x;  As=Ax; ATAs=ATAx;
            end

            g=ATAs-ATy + lambda * ones(n,1);
            % compute the gradient g

            v=s-g/L;
            % a gradient step based on the search point s

            % projection, we need to revise this part
            switch(opts.treeFlag)
                case 1  
                    x=sequence_bottomup(v,n);
                case 2  
                    x=orderTreeBinary(v,n);
                case 3  
                    x=orderTreeDepth1(v,n);
                case 4  
                    x=orderTree(FileName, v, rootNum, n);
                otherwise
                    x=sequence_bottomup(v,n);
            end

            v=xnew-s;  % the difference between the new approximate solution x
            % and the search point s

            % compute A xnew
            if (opts.nFlag==0)
                Axnew=A* xnew;
            elseif (opts.nFlag==1)
                invNu=xnew./nu; mu_invNu=mu * invNu;
                Axnew=A*invNu -repmat(mu_invNu, m, 1);
            else
                Axnew=A*xnew-repmat(mu*xnew, m, 1);     Axnew=Axnew./nu;
            end

            Av=Axnew -As;
            r_sum=v'*v; l_sum=Av'*Av;
            
            if (r_sum <=1e-20)
                bFlag=1; % this shows that, the gradient step makes little improvement
                break;
            end

            % the condition is ||Av||_2^2 <= L* ||v||_2^2
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

        Axy=Ax-y;
        funVal(iterStep)=Axy' * Axy/2 + lambda* sum(x);
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
end

