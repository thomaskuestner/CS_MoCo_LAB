function res=tesla(A, y, lambda1, lambda2, opts)
%
%%
% Function mtLogisticR:
%      Logistic Loss for ...... (to be added)
%
%% Problem
%
%   to be added
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size m x 1)
%  lambda1 -
%  lambda2 -
%  opts-      Optional inputs (default value: opts=[])
%
%% Output parameters:
%  x-         The obtained weight of size n x k
%  c-         The obtained intercept if size 1 x k
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
% Last modified 3 November 2009.
%
% Related functions:
%  ....
%%

%% Verify and initialize the parameters
%%
if (nargin <5)
    error('\n Inputs: A, y, lambda1,  lambda2, and opts.ind should be specified!\n');
end

[m,n]=size(A);

if (length(y) ~=m)
    error('\n Check the length of y!\n');
end

if (lambda1<0 || lambda2 <0)
    error('\n lambda1 should be nonnegative!\n');
end

opts=sll_opts(opts); % run sll_opts to set default values (flags)

%% Detailed initialization
%%

% Initialize ind and q
if ~isfield(opts,'ind')
    error('\n In tesla, .ind should be specified');
else
    ind=opts.ind;
    k=length(ind)-1;

    if ind(k+1)~=m
        error('\n Check opts.ind');
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

p_flag=(y==1);                  % the indices of the postive samples

for i=1:k
    ind_i=(ind(i)+1):ind(i+1);  % indices for the i-th group

    m1(1,i)=sum(p_flag(ind_i));      % the total number of the positive samples
    m2(1,i)=length(ind_i)-m1(1,i);   % the total number of the negative samples
end

% process the regularization parameter
if (opts.rFlag==1) % lambda1 here is the scaling factor lying in [0,1]
    if (lambda1<0 || lambda1>1)
        error('\n opts.rFlag=1, and lambda1 should be in [0,1]');
    end

    % we compute ATb for computing lambda_max, when the input lambda1 is a ratio

    p_flag=(y==1);                  % the indices of the postive samples

    for i=1:k
        ind_i=(ind(i)+1):ind(i+1);  % indices for the i-th group

        b(ind_i,1)=p_flag(ind_i)*m1(i)/(m1(i)+m2(i))-...
            (~p_flag(ind_i))*m2(i)/(m1(i)+m2(i));
    end

    ATb=zeros(n, k);
    % compute AT b
    for i=1:k
        ind_i=(ind(i)+1):ind(i+1);     % indices for the i-th group

        if (opts.nFlag==0)
            tt =A(ind_i,:)'*b(ind_i,1);
        elseif (opts.nFlag==1)
            tt= A(ind_i,:)'*b(ind_i,1) - sum(b(ind_i,1)) * mu';
            tt=tt./nu(ind_i,1);
        else
            invNu=b(ind_i,1)./nu(ind_i,1);
            tt=A(ind_i,:)'*invNu - sum(invNu)*mu';
        end

        ATb(:,i)= tt;
    end

    q_bar=Inf;
    
    lambda_max=0;
    for i=1:n
        lambda_max=max(lambda_max,...
            norm(  ATb(i,:), q_bar) );
    end

    lambda_max=lambda_max / m;
    lambda1=lambda1*lambda_max;
    
    if isfield(opts,'rFlag2')
        if (opts.rFlag2==1)
            lambda2=lambda1*lambda_max;
        end
    end
end

% initialize a starting point
if opts.init==2
    x=zeros(n,k); c=zeros(1,k);
else
    if isfield(opts,'x0')
        x=opts.x0;
        if ( size(x,1)~=n && size(x,2)~=k )
            error('\n Check the input .x0');
        end
    else
        x=zeros(n,k);
    end

    if isfield(opts,'c0')
        c=opts.c0;

        if ( length(c)~=k )
            error('\n Check the input .c0');
        end
    else
        c=log(m1./m2);
    end
end

Ax=zeros(m,1); % m x 1
% compute Ax: Ax_i= A_i * x_i
for i=1:k
    ind_i=(ind(i)+1):ind(i+1);     % indices for the i-th group
    m_i=ind(i+1)-ind(i);          % number of samples in the i-th group

    if (opts.nFlag==0)
        Ax(ind_i,1)=A(ind_i,:)* x(:,i);
    elseif (opts.nFlag==1)
        invNu=x(:,i)./nu; mu_invNu=mu * invNu;
        Ax(ind_i,1)=A(ind_i,:)*invNu -repmat(mu_invNu, m_i, 1);
    else
        Ax(ind_i,1)=A(ind_i,:)*x(:,i)-repmat(mu*x(:,i), m, 1);
        Ax(ind_i,1)=Ax./nu(ind_i,1);
    end
end


%% The main program
% The Armijo Goldstein line search schemes + accelearted gradient descent

z0=zeros(n-1,k);

if (opts.lFlag==0)

    bFlag=0; % this flag tests whether the gradient step only changes a little

    L=1/m; % the intial guess of the Lipschitz continuous gradient

    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,k);
    cp=c;         ccp=zeros(1,k);

    alphap=0; alpha=1;

    for iterStep=1:opts.maxIter
        % --------------------------- step 1 ---------------------------
        % compute search point s based on xp and x (with beta)
        beta=(alphap-1)/alpha;    s=x + beta* xxp;   sc=c + beta* ccp;

        % --------------------------- step 2 ---------------------------
        % line search for L and compute the new approximate solution x

        % compute the gradient (g) at s
        As=Ax + beta* (Ax-Axp);

        % aa= - diag(y) * (A * s + sc)
        vec_sc=zeros(m,1);
        for i=1:k
            ind_i=(ind(i)+1):ind(i+1);     % indices for the i-th group
            vec_sc(ind_i,1)=sc(i);
        end
        aa=- y.*(As+ vec_sc);

        % fun_s is the logistic loss at the search point
        bb=max(aa,0);
        fun_s= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / m;

        % compute prob=[p_1;p_2;...;p_m]
        prob=1./( 1+ exp(aa) );

        % b= - diag(y) * (1 - prob)
        b= -y.*(1-prob) / m;

        gc=zeros(1,k);                     % the gradient of c
        for i=1:k
            ind_i=(ind(i)+1):ind(i+1);     % indices for the i-th group
            gc(1,k)=sum(b(ind_i));
        end

        % compute g= AT b, the gradient of x
        for i=1:k
            ind_i=(ind(i)+1):ind(i+1);     % indices for the i-th group

            if (opts.nFlag==0)
                tt =A(ind_i,:)'*b(ind_i,1);
            elseif (opts.nFlag==1)
                tt= A(ind_i,:)'*b(ind_i,1) - sum(b(ind_i,1)) * mu';
                tt=tt./nu(ind_i,1);
            else
                invNu=b(ind_i,1)./nu(ind_i,1);
                tt=A(ind_i,:)'*invNu - sum(invNu)*mu';
            end

            g(:,i)= tt;
        end

        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;
        cp=c;

        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and then do the Lq/L1-norm regularized projection
            v=s-g/L; c= sc- gc/L;

            [x, z, gap]=tesla_proj(v, z0, ...
                lambda1/L, lambda2/L, n, k,...
                1000, 1e-8, 1, 6);
            z0=z;
            
            v=x-s;  % the difference between the new approximate solution x
            % and the search point s

            % compute Ax: Ax_i= A_i * x_i
            for i=1:k
                ind_i=(ind(i)+1):ind(i+1);     % indices for the i-th group
                m_i=ind(i+1)-ind(i);          % number of samples in the i-th group

                if (opts.nFlag==0)
                    Ax(ind_i,1)=A(ind_i,:)* x(:,i);
                elseif (opts.nFlag==1)
                    invNu=x(:,i)./nu; mu_invNu=mu * invNu;
                    Ax(ind_i,1)=A(ind_i,:)*invNu -repmat(mu_invNu, m_i, 1);
                else
                    Ax(ind_i,1)=A(ind_i,:)*x(:,i)-repmat(mu*x(:,i), m, 1);
                    Ax(ind_i,1)=Ax./nu(ind_i,1);
                end
            end

            % aa= - diag(y) * (A * x + c)
            vec_sc=zeros(m,1);
            for i=1:k
                ind_i=(ind(i)+1):ind(i+1);     % indices for the i-th group
                vec_sc(ind_i,1)=c(i);
            end
            aa=- y.*(Ax+ vec_sc);

            % fun_s is the logistic loss at the search point
            bb=max(aa,0);
            fun_x= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / m;

            r_sum=(norm(v,'fro')^2 + norm(c-sc,2)^2) / 2;
            l_sum=fun_x - fun_s - sum(sum(v.* g)) - (c-sc)* gc';

            if (r_sum <=1e-20)
                bFlag=1; % this shows that, the gradient step makes little improvement
                break;
            end

            % the condition is fun_x <= fun_s + <v, g> + <c ,gc>
            %                           + L/2 * (<v,v> + <c-sc,c-sc> )
            if(l_sum <= r_sum * L)
                break;
            else
                L=max(2*L, l_sum/r_sum);
                % fprintf('\n L=%5.6f',L);
            end
        end

        % --------------------------- step 3 ---------------------------
        % update alpha and alphap, and check whether converge
        alphap=alpha; alpha= (1+ sqrt(4*alpha*alpha +1))/2;

        ValueL(iterStep)=L;

        xxp=x-xp;   ccp=c-cp;

        funVal(iterStep)=fun_x;

        funVal(iterStep)=funVal(iterStep)+ lambda1 * sum(abs(x(:))) + lambda2 * sum( sum( abs(x(:,1:(k-1)) - x(:,2:k) ) ) );

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
        res.gap(:,iterStep)=gap;
    end
    res.x=x;
    res.c=c;
    res.funVal=funVal;
    res.ValueL=ValueL;
    res.lambda_max=lambda_max;
end
