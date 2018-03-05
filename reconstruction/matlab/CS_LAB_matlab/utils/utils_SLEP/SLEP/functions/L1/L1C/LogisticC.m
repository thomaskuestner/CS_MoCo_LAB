function [x, c, funVal, ValueL]=LogisticC(A, y, z, opts)
%
%%
% Function LogisticC
%      Logistic Loss with the L1 ball constraint
%
%% Problem
%
%  min  f(x,c) = \sum_i - weight_i * log (p_i) + 1/2 rsL2 * ||x||_2^2
%  s.alpha. \|x\|_1 <=z
%
%  a_i denotes a training sample,
%      and a_i' corresponds to the i-th row of the data matrix A
%
%  y_i (either 1 or -1) is the response
%
%  p_i= 1/ (1+ exp(-y_i (x' * a_i + c) ) ) denotes the probability
%
%  weight_i denotes the weight for the i-th sample
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size mx1)
%  z -        Radius of the L1 ball
%  opts-      Optional inputs (default value: opts=[])
%
%% Output parameters:
%  x-         The obtained weight of size n x 1
%  c-         The obtained intercept (scalar)
%  funVal-    Function value during iterations
%
%% Copyright (C) 2009-2010 Jun Liu, and Jieping Ye
%
% You are suggested to first read the Manual.
%
% For any problem, please contact with Jun Liu via j.liu@asu.edu
%
%% Related papers
%
% [1]  Jun Liu, Jianhui Chen, and Jieping Ye, Large-Scale Sparse Logistic
%      Regression, KDD 2009
%
% [2]  Jun Liu and Jieping Ye, Efficient Euclidean Projections
%      in Linear Time, ICML 2009.
%
%% Related functions
%
%  sll_opts, initFactor, pathSolutionLeast,
%  LogisticR, nnLogisticR, nnLogisticC,
%  eplb
%
%%

%% Verify and initialize the parameters
%%
if (nargin <3)
    error('\n Inputs: A, y and z should be specified!\n');
end

[m,n]=size(A);

if (length(y) ~=m)
    error('\n Check the length of y!\n');
end

if (z<=0)
    error('\n z should be positive!\n');
end

opts=sll_opts(opts); % run sll_opts to set default values (flags)

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

%% Group & Others

% The parameter 'weight' contains the weight for each training sample.
% See the definition of the problem above.
% The summation of the weights for all the samples equals to 1.
if (isfield(opts,'sWeight'))
    sWeight=opts.sWeight;

    if ( length(sWeight)~=2 || sWeight(1) <=0 || sWeight(2) <= 0)
        error('\n Check opts.sWeight, which contains two positive values');
    end

    % we process the weight, so that the summation of the weights for all
    % the samples is 1.

    p_flag=(y==1);                  % the indices of the postive samples
    m1=sum(p_flag) * sWeight(1);    % the total weight for the positive samples
    m2=sum(~p_flag) * sWeight(2);   % the total weight for the positive samples

    weight(p_flag,1)=sWeight(1)/(m1+m2);
    weight(~p_flag,1)=sWeight(2)/(m1+m2);
else
    weight=ones(m,1)/m;             % if not specified, we apply equal weight
end

%% Starting point initialization

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
    x=zeros(n,1); c=0;
else
    if isfield(opts,'x0')
        x=opts.x0;
        if (length(x)~=n)
            error('\n Check the input .x0');
        end
    else
        x=zeros(n,1);
    end

    if isfield(opts,'c0')
        c=opts.c0;
    else
        c=0;
    end
end

% compute A x
if (opts.nFlag==0)
    Ax=A* x;
elseif (opts.nFlag==1)
    invNu=x./nu; mu_invNu=mu * invNu;
    Ax=A*invNu -repmat(mu_invNu, m, 1);
else
    Ax=A*x-repmat(mu*x, m, 1);    Ax=Ax./nu;
end

%% The main program

%% Nemirovski's line search

if (opts.lFlag==0)

    lambda0=0;
    % a guess of the root

    bFlag=0; % this flag tests whether the gradient step only changes a little

    L=1/m + rsL2; % the intial guess of the Lipschitz continuous gradient

    weighty=weight.*y;
    % the product between weight and y

    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,1);
    cp=c;         ccp=0;

    %% The Armijo Goldstein line search schemes + accelearted gradient descent

    alphap=0; alpha=1;

    for iterStep=1:opts.maxIter
        % --------------------------- step 1 ---------------------------
        % compute search point s based on xp and x (with beta)
        beta=(alphap-1)/alpha;    s=x + beta* xxp;  sc=c + beta* ccp;

        % --------------------------- step 2 ---------------------------
        % line search for L and compute the new approximate solution x

        % compute As=A*s
        As=Ax + beta* (Ax-Axp);

        % aa= - diag(y) * (A * s + sc)
        aa=- y.*(As+ sc);

        % fun_s is the logistic loss at the search point
        bb=max(aa,0);
        fun_s= weight' * ( log( exp(-bb) +  exp(aa-bb) ) + bb )+...
            rsL2/2 * s'*s;

        % compute prob=[p_1;p_2;...;p_m]
        prob=1./( 1+ exp(aa) );

        % b= - diag(y.* weight) * (1 - prob)
        b= -weighty.*(1-prob);

        gc=sum(b); % the gradient of c

        % compute g= AT b, the gradient of x
        if (opts.nFlag==0)
            g=A'*b;
        elseif (opts.nFlag==1)
            g=A'*b - sum(b) * mu';  g=g./nu;
        else
            invNu=b./nu;              g=A'*invNu-sum(invNu)*mu';
        end

        g=g+ rsL2 * s;  % add the squared L2 norm regularization term

        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;
        cp=c;

        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and then do the Lq/L1-norm regularized projection
            v=s-g/L;  c= sc- gc/L;

            % projection
            [x, lambda, zf_step]=eplb(v, n, z, lambda0);
            lambda0=lambda;

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

            % aa= - diag(y) * (A * x + c)
            aa=- y.*(Ax+ c);

            % fun_x is the logistic loss at the new approximate solution
            bb=max(aa,0);
            fun_x= weight'* ( log( exp(-bb) +  exp(aa-bb) ) + bb )+...
                rsL2/2 * x'*x;

            r_sum= (v'*v + (c-sc)^2) / 2;
            l_sum=fun_x - fun_s - v'* g - (c-sc)* gc;

            if (r_sum <=1e-20)
                bFlag=1; % this shows that, the gradient step makes little improvement
                break;
            end

            % the condition is fun_x <= fun_s + v'* g + c * gc
            %                           + L/2 * (v'*v + (c-sc)^2 )
            if(l_sum <= r_sum * L)
                break;
            else
                L=max(2*L, l_sum/r_sum);
                %fprintf('\n L=%e, r_sum=%e',L, r_sum);
            end
        end

        % --------------------------- step 3 ---------------------------
        % update alpha and alphap, and check whether converge
        alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;
        
        ValueL(iterStep)=L;
        % store values for L
        
        xxp=x-xp;    ccp=c-cp;
        funVal(iterStep)=fun_x;

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

%% adaptive line search

if (opts.lFlag==1)

    lambda0=0;
    % a guess of the root

    bFlag=0; % this flag tests whether the gradient step only changes a little

    L=1/m + rsL2; % the intial guess of the Lipschitz continuous gradient
    
    gamma=1;
    % we shall set the value of gamma = L,
    % and L is appropriate for the starting point

    weighty=weight.*y;
    % the product between weight and y

    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,1);
    cp=c;         ccp=0;

    for iterStep=1:opts.maxIter
        while (1)
            if (iterStep~=1)
                alpha= (-gamma+ sqrt(gamma*gamma + 4* L * gamma)) / (2*L);
                beta= (gamma - gamma* alphap) / (alphap * gamma + alphap* L * alpha);
                % beta is the coefficient for generating search point s

                s=x + beta* xxp;  sc=c + beta* ccp;
                As=Ax + beta* (Ax-Axp);
            else
                alpha= (-1+ sqrt(5)) / 2;  beta=0;
                s=x; sc= c;
                As=Ax;
            end

            % aa= - diag(y) * (A * s + sc)
            aa=- y.*(As+ sc);

            % fun_s is the logistic loss at the search point
            bb=max(aa,0);
            fun_s= weight' * ( log( exp(-bb) +  exp(aa-bb) ) + bb )+...
                rsL2/2 * s'*s;

            % compute prob=[p_1;p_2;...;p_m]
            prob=1./( 1+ exp(aa) );

            % b= - diag(y.* weight) * (1 - prob)
            b= -weighty.*(1-prob);

            gc=sum(b); % the gradient of c

            % compute g= AT b, the gradient of x
            if (opts.nFlag==0)
                g=A'*b;
            elseif (opts.nFlag==1)
                g=A'*b - sum(b) * mu';  g=g./nu;
            else
                invNu=b./nu;              g=A'*invNu-sum(invNu)*mu';
            end

            g=g+ rsL2 * s;  % add the squared L2 norm regularization term

            % let s walk in a step in the antigradient of s to get v
            % and then do the Lq/L1-norm regularized projection
            v=s-g/L;  cnew= sc- gc/L;

            % projection
            [xnew, lambda, zf_step]=eplb(v, n, z, lambda0);
            lambda0=lambda;

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

            % aa= - diag(y) * (A * x + c)
            aa=- y.*(Axnew+ cnew);

            % fun_x is the logistic loss at the new approximate solution
            bb=max(aa,0);
            fun_x= weight'* ( log( exp(-bb) +  exp(aa-bb) ) + bb )+...
                rsL2/2 * xnew'*xnew;

            r_sum= (v'*v + (cnew-sc)^2) / 2;
            l_sum=fun_x - fun_s - v'* g - (cnew-sc)* gc;

            if (r_sum <=1e-20)
                bFlag=1; % this shows that, the gradient step makes little improvement
                break;
            end

            % the condition is fun_x <= fun_s + v'* g + c * gc
            %                           + L/2 * (v'*v + (c-sc)^2 )
            if(l_sum <= r_sum * L)
                break;
            else
                L=max(2*L, l_sum/r_sum);
                %fprintf('\n L=%e, r_sum=%e',L, r_sum);
            end
        end

        gamma=L* alpha* alpha;    alphap=alpha;
        % update gamma, and alphap

        ValueL(iterStep)=L;
        % store values for L

        tao=L * r_sum / l_sum;
        if (tao >=5)
            L=L*0.8;
        end
        % decrease the value of L
        
        xp=x;     x=xnew; xxp=x-xp;
        Axp=Ax;   Ax=Axnew;                
        cp=c;     c=cnew; ccp=c-cp;
        % update x and Ax with xnew and Axnew
        
        funVal(iterStep)=fun_x;

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