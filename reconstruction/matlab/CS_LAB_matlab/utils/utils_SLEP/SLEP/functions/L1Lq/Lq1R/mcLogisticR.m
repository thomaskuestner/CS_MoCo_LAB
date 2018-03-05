function [x, c, funVal, ValueL]=mcLogisticR(A, y, z, opts)
%
%%
% Function mcLogisticR:
%      Logistic Loss for Multi-class (task) Learning
%             via the (group) L1/Lq-norm Regularization
%
%% Problem
%
%  min  - sum_{il} weight_{il} log( p_{il} ) + z * sum_j ||x^j||_q
%
%  p_{il}= 1 / (1+ exp(-y_i (x_i' * a_i + c_l) ) ) denotes the probability
%  weight_{il} is the weight for the i-th sample in the l-th classifier
%                is a m x k matrix
%  c_l is the intercept for the l-th classfier, and is a 1xk vector
%  x_i denotes the i-th column of x
%  x^j denotes the j-th row of x
%  a_i' denotes the i-th row of A
%
%  y_i (either 1 or -1) is the response
%
%  In this implementation, we assume weight_{il}=1/(mk)
%
%  When y are of binary values (1 and -1), this problem is the
%  multi-class classificaiton problem, which learns the projection
%  matrix x, by sharing information across k classification tasks (with the
%  L1/Lq-norm Regularization)
%
%  For the case that the multi tasks have different data
%  matrix, please refer to the functions:
%            mtLeastR and mtLogisticR.
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size m x k)
%  z -        L1/Lq norm regularization parameter (z >=0)
%  opts-      Optional inputs (default value: opts=[])
%
%% Output parameters:
%  x-         The obtained weight of size n x k
%  c-         The obtained intercept if size 1 x k
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
%  mcLeastR, eppMatrix
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


%% Group Property

% Initialize q
if (~isfield(opts,'q'))
    q=2; opts.q=2;
else
    q=opts.q;
    if (q<1)
        error('\n q should be larger than 1');
    end
end


%% Starting point initialization

% % % compute AT y
% % if (opts.nFlag==0)
% %     ATy =A'* y;
% % elseif (opts.nFlag==1)
% %     ATy= A'* y -  mu' * sum(y, 1);    ATy=ATy./repmat(nu, 1, k);
% % else
% %     invNu=y./repmat(nu, 1, k);        ATy=A'*invNu-mu' * sum(invNu, 1);
% % end


p_flag=(y==1);                  % the indices of the postive samples
m1=sum(p_flag,1);               % the total number of the positive samples
m2=m-m1;                        % the total number of the negative samples
% m1 and m2 are 1 x k vectors

% process the regularization parameter
if (opts.rFlag==0)
    lambda=z;
else % z here is the scaling factor lying in [0,1]
    if (z<0 || z>1)
        error('\n opts.rFlag=1, and z should be in [0,1]');
    end
    
    % we compute ATb for computing lambda_max, when the input z is a ratio
    
    b=zeros(m, k);                  % b is of size m x k
    b= b + p_flag.*repmat(m2,m,1);
    b= b -(~p_flag).*repmat(m1,m,1);
    b=b/m;
    
    % compute AT b
    if (opts.nFlag==0)
        ATb =A'* b;
    elseif (opts.nFlag==1)
        ATb= A'* b -  mu' * sum(b, 1);    ATb=ATb./repmat(nu, 1, k);
    else
        invNu=b./repmat(nu, 1, k);        ATb=A'*invNu-mu' * sum(invNu, 1);
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
            norm(  ATb(i,:), q_bar) );
    end
    
    lambda_max=lambda_max / (m*k);
    lambda=z*lambda_max;
end

% initialize a starting point
if opts.init==2
    x=zeros(n,k); c=log(m1./m2);
else
    if isfield(opts,'x0')
        x=opts.x0;
        if ( size(x,1)~=n || size(x,2)~=k )
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

% compute Ax= A * x
if (opts.nFlag==0)
    Ax=A* x;
elseif (opts.nFlag==1)
    invNu=x./repmat(nu, 1, k);  mu_invNu=mu * invNu;
    Ax=A*invNu -repmat(mu_invNu, m, 1);
else
    Ax=A*x-repmat(mu*x, m, 1);     Ax=Ax./repmat(nu, 1, k);
end

%% The main program

if (opts.mFlag==0 && opts.lFlag==0)
    % The Armijo Goldstein line search schemes + accelearted gradient descent
    
    bFlag=0; % this flag tests whether the gradient step only changes a little
    
    L=1/(m*k); % the intial guess of the Lipschitz continuous gradient
    
    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,k);
    cp=c;         ccp=zeros(1,k);
    
    alphap=0; alpha=1;
    
    for iterStep=1:opts.maxIter
        % --------------------------- step 1 ---------------------------
        % compute search point s based on xp and x (with beta)
        beta=(alphap-1)/alpha;    s=x + beta* xxp;  sc=c + beta* ccp;
        
        % --------------------------- step 2 ---------------------------
        % line search for L and compute the new approximate solution x
        
        % compute the gradient (g) at s
        As=Ax + beta* (Ax-Axp);
        
        % aa= - diag(y) * (A * s + sc)
        aa=- y.*(As+ repmat(sc,m,1));
        
        % fun_s is the logistic loss at the search point
        bb=max(aa,0);
        fun_s= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / (m * k);
        
        % compute prob=[p_1;p_2;...;p_m]
        prob=1./( 1+ exp(aa) );
        
        % b= - diag(y) * (1 - prob)
        b= -y.*(1-prob) / (m * k);
        
        gc=sum(b); % the gradient of c
        
        % compute g= AT b, the gradient of x
        if (opts.nFlag==0)
            g=A'*b;
        elseif (opts.nFlag==1)
            g= A'* b -  mu' * sum(b, 1);  g=g./repmat(nu, 1, k);
        else
            invNu=b./repmat(nu, 1, k);        g=A'*invNu-mu' * sum(invNu, 1);
        end
        
        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;
        cp=c;
        
        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and then do the L1/Lq-norm regularized projection
            v=s-g/L; c= sc- gc/L;
            
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
            
            % aa= - diag(y) * (A * x + c)
            aa=- y.*(Ax+ repmat(c,m,1));
            
            % fun_x is the logistic loss at the new approximate solution
            bb=max(aa,0);
            fun_x= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / (m * k);
            
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
        
        xxp=x-xp;  ccp=c-cp;
        
        funVal(iterStep)=fun_x;
        
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
    end
end


%%
% The Reformulated problem + Armijo Goldstein line search scheme
%
% Original Problem:
%  min  - sum_{il} weight_{il} log( p_{il} ) + z * sum_j ||x^j||_2
%
% Reformulated Problem:
%  min  - sum_{il} weight_{il} log( p_{il} ) + z * 1^ t
%  s.t.   ||x^j||_2 <= t_j


if (opts.mFlag==1 && opts.lFlag==0 && opts.q==2)
    
    bFlag=0; % this flag tests whether the gradient step only changes a little
    
    L=1/(m*k); % the intial guess of the Lipschitz continuous gradient
    
    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,k);
    cp=c;         ccp=zeros(1,k);
    
    t=sqrt(sum(x.^2, 2)); tp=t;
    % t is the upper bound of absolute value of x
    
    alphap=0; alpha=1;
    
    for iterStep=1:opts.maxIter
        % --------------------------- step 1 ---------------------------
        % compute search point s based on xp and x (with beta)
        beta=(alphap-1)/alpha;    s=x + beta* xxp;  sc=c + beta* ccp;
        s_t= t + beta * (t-tp);
        
        % --------------------------- step 2 ---------------------------
        % line search for L and compute the new approximate solution x
        
        % compute the gradient (g) at s
        As=Ax + beta* (Ax-Axp);
        
        % aa= - diag(y) * (A * s + sc)
        aa=- y.*(As+ repmat(sc,m,1));
        
        % fun_s is the logistic loss at the search point
        bb=max(aa,0);
        fun_s= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / (m * k);
        
        % compute prob=[p_1;p_2;...;p_m]
        prob=1./( 1+ exp(aa) );
        
        % b= - diag(y) * (1 - prob)
        b= -y.*(1-prob) / (m * k);
        
        gc=sum(b); % the gradient of c
        
        % compute g= AT b, the gradient of x
        if (opts.nFlag==0)
            g=A'*b;
        elseif (opts.nFlag==1)
            g= A'* b -  mu' * sum(b, 1);  g=g./repmat(nu, 1, k);
        else
            invNu=b./repmat(nu, 1, k);        g=A'*invNu-mu' * sum(invNu, 1);
        end
        
        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;
        cp=c;    tp=t;
        
        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and then do the L1/Lq-norm regularized projection
            u=s-g/L; c= sc- gc/L; v= s_t - lambda / L;
            
            % L1/Lq-norm regularized projection
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
            
            % aa= - diag(y) * (A * x + c)
            aa=- y.*(Ax+ repmat(c,m,1));
            
            % fun_x is the logistic loss at the new approximate solution
            bb=max(aa,0);
            fun_x= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / (m * k);
            
            r_sum=(norm(v,'fro')^2 + norm(c-sc,2)^2 + norm(t-s_t,2)^2) / 2;
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
        
        xxp=x-xp;  ccp=c-cp;
        
        funVal(iterStep)=fun_x + lambda * sum(t);
        
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
% The Reformulated problem + Adaptive line search 
%
% Original Problem:
%  min  - sum_{il} weight_{il} log( p_{il} ) + z * sum_j ||x^j||_2
%
% Reformulated Problem:
%  min  - sum_{il} weight_{il} log( p_{il} ) + z * 1^ t
%  s.t.   ||x^j||_2 <= t_j


if (opts.mFlag==1 && opts.lFlag==1 && opts.q==2)
    
    bFlag=0; % this flag tests whether the gradient step only changes a little
    
    L=1/(m*k); % the intial guess of the Lipschitz continuous gradient
    
    gamma=1;
    % we shall set the value of gamma = L,
    % and L is appropriate for the starting point
    
    % assign xp with x, and Axp with Ax
    xp=x; Axp=Ax; xxp=zeros(n,k);
    cp=c;         ccp=zeros(1,k);
    t=sqrt(sum(x.^2, 2)); tp=t;
    % t is the upper bound of absolute value of x
    
    for iterStep=1:opts.maxIter
        while (1)            
            if (iterStep~=1)
                alpha= (-gamma+ sqrt(gamma*gamma + 4* L * gamma)) / (2*L);
                beta= (gamma - gamma* alphap) / (alphap * gamma + alphap* L * alpha);
                % beta is the coefficient for generating search point s
                
                s=x + beta* xxp;   s_t= t + beta * (t -tp); sc=c + beta* ccp;
                As=Ax + beta* (Ax-Axp);
            else
                alpha= (-1+ sqrt(5)) / 2;  beta=0; 
                s=x; s_t=t; sc= c; 
                As=Ax;
            end
            
            % aa= - diag(y) * (A * s + sc)
            aa=- y.*(As+ repmat(sc,m,1));
            
            % fun_s is the logistic loss at the search point
            bb=max(aa,0);
            fun_s= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / (m * k);
            
            % compute prob=[p_1;p_2;...;p_m]
            prob=1./( 1+ exp(aa) );
            
            % b= - diag(y) * (1 - prob)
            b= -y.*(1-prob) / (m * k);
            
            gc=sum(b); % the gradient of c
            
            % compute g= AT b, the gradient of x
            if (opts.nFlag==0)
                g=A'*b;
            elseif (opts.nFlag==1)
                g= A'* b -  mu' * sum(b, 1);  g=g./repmat(nu, 1, k);
            else
                invNu=b./repmat(nu, 1, k);        g=A'*invNu-mu' * sum(invNu, 1);
            end
            
            % let s walk in a step in the antigradient of s to get v
            % and then do the L1/Lq-norm regularized projection
            u=s-g/L; cnew= sc- gc/L; v= s_t - lambda / L;
            
            % L1/Lq-norm regularized projection
            [xnew, tnew]=ep21R(u,v,n,k);
            
            v=xnew-s;  % the difference between the new approximate solution x
            % and the search point s
            
            % compute A * xnew
            if (opts.nFlag==0)
                Axnew=A* xnew;
            elseif (opts.nFlag==1)
                invNu=xnew./nu; mu_invNu=mu * invNu;
                Axnew=A*invNu -repmat(mu_invNu, m, 1);
            else
                Axnew=A*xnew-repmat(mu*xnew, m, 1);     Axnew=Axnew./nu;
            end
            
            % aa= - diag(y) * (A * xnew + cnew)
            aa=- y.*(Axnew+ repmat(cnew,m,1));
            
            % fun_x is the logistic loss at the new approximate solution
            bb=max(aa,0);
            fun_x= sum(sum ( log( exp(-bb) +  exp(aa-bb) ) + bb ) ) / (m * k);
            
            r_sum=(norm(v,'fro')^2 + norm(cnew-sc,2)^2 + norm(tnew-s_t,2)^2) / 2;
            l_sum=fun_x - fun_s - sum(sum(v.* g)) - (cnew-sc)* gc';
            
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
        
        gamma=L* alpha* alpha;    alphap=alpha;
        % update gamma, and alphap

        ValueL(iterStep)=L;
        % store values for L

        tao=L * r_sum / l_sum;
        if (tao >=5)
            L=L*0.8;
        end
        % decrease the value of L

        xp=x;     Axp=Ax;
        x=xnew;   Ax=Axnew;
        tp=t; t=tnew;        
        cp=c; c=cnew;
        % update x, t, c, and Ax        
        xxp=x-xp;  ccp=c-cp;
        
        funVal(iterStep)=fun_x + lambda * sum(t);
        
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