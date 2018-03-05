function [x, c, funVal, ValueL]=tree_mcLogisticR(A, y, z, opts)
%
%%
% Function mcLogisticR:
%      Logistic Loss for Multi-class (task) Learning
%             with the tree structured group Lasso Regularization
%
%% Problem
%
%  min  - sum_{il} weight_{il} log( p_{il} ) + z * sum_i sum_j w_j ||x^i_{G_j}||
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
%  x^i denotes the i-th row of x
%
%  G_j's are nodes with tree structure
%
%  The tree structured group information is contained in
%  opts.ind, which is a 3 x nodes matrix, where nodes denotes the number of
%  nodes of the tree.
%
%  opts.ind(1,:) contains the starting index
%  opts.ind(2,:) contains the ending index
%  opts.ind(3,:) contains the corresponding weight (w_j)
%
%  Note: 
%  1) If each element of x^j is a leaf node of the tree and the weight for
%  this leaf node are the same, we provide an alternative "efficient" input
%  for this kind of node, by creating a "super node" with 
%  opts.ind(1,1)=-1; opts.ind(2,1)=-1; and opts.ind(3,1)=the common weight.
%
%  2) If the features are well ordered in that, the features of the left
%  tree is always less than those of the right tree, opts.ind(1,:) and
%  opts.ind(2,:) contain the "real" starting and ending indices. That is to
%  say, x^j( opts.ind(1,j):opts.ind(2,j) ) denotes x^j_{G_j}. In this case,
%  the entries in opts.ind(1:2,:) are within 1 and k (the number of tasks).
%
%
%  If the features are not well ordered, please use the input opts.G for
%  specifying the index so that  
%   x^j( opts.G ( opts.ind(1,j):opts.ind(2,j) ) ) denotes x^j_{G_j}.
%  In this case, the entries of opts.G are within 1 and n, and the entries of
%  opts.ind(1:2,:) are within 1 and length(opts.G).
%
% The following example shows how G and ind works:
%
% G={ {1, 2}, {4, 5}, {3, 6}, {7, 8},
%     {1, 2, 3, 6}, {4, 5, 7, 8}, 
%     {1, 2, 3, 4, 5, 6, 7, 8} }.
%
% ind={ [1, 2, 100]', [3, 4, 100]', [5, 6, 100]', [7, 8, 100]',
%       [9, 12, 100]', [13, 16, 100]', [17, 24, 100]' },
%
% where each node has a weight of 100.
%
%% Input parameters:
%
%  A-         Matrix of size m x n
%                A can be a dense matrix
%                         a sparse matrix
%                         or a DCT matrix
%  y -        Response vector (of size m x k)
%  z -        Regularization parameter (z >=0)
%  opts-      optional inputs (default value: opts=[])
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
% Last modified on October 3, 2010.
%
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
if (nargin <4)
    error('\n Inputs: A, y, z, and opts (.ind) should be specified!\n');
end

[m,n]=size(A);

k=size(y,2);

if (size(y,1) ~=m)
    error('\n Check the length of y!\n');
end

if (z<0)
    error('\n z should be nonnegative!\n');
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


%% Group and others

% Initialize ind 
if (~isfield(opts,'ind'))
    error('\n In tree_mcLeastR, the field .ind should be specified');
else
    ind=opts.ind;
   
    if (size(ind,1)~=3)
        error('\n Check opts.ind');
    end
end


GFlag=0;
% if GFlag=1, we shall apply general_altra
if (isfield(opts,'G'))
    GFlag=1;
    
    G=opts.G;
    if (max(G) >k || max(G) <1)
        error('\n The input G is incorrect. It should be within %d and %d',1,k);
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
%     if (z<0 || z>1)
%         error('\n opts.rFlag=1, and z should be in [0,1]');
%     end
    
    % we compute ATb for computing lambda_max, when the input z is a ratio
    
    computedFlag=0;
    if (isfield(opts,'lambdaMax'))
        if (opts.lambdaMax~=-1)
            lambda=z*opts.lambdaMax;
            computedFlag=1;
        end
    end

    if (~computedFlag)
        
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

        if (GFlag==0)
            lambda_max=findLambdaMax_mt(ATb, n, k, ind, size(ind,2));
        else
            lambda_max=general_findLambdaMax_mt(ATb, n, k, G, ind, size(ind,2));
        end
        
        lambda_max= lambda_max / (m*k);
        
        % As .rFlag=1, we set lambda as a ratio of lambda_max
        lambda=z*lambda_max;
    end
end

% The following is for computing lambdaMax
% we use opts.lambdaMax=-1 to show that we need the computation.
% 
% One can use this for setting up opts.lambdaMax
if (isfield(opts,'lambdaMax'))
    
    if (opts.lambdaMax==-1)
        
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
        
        % compute lambda_max
        if (GFlag==0)
            lambda_max=findLambdaMax_mt(ATb, n, k, ind, size(ind,2));
        else
            lambda_max=general_findLambdaMax_mt(ATb, n, k, G, ind, size(ind,2));
        end
        
        x=lambda_max;
        funVal=lambda_max;
        ValueL=lambda_max;
        
        return;
    end
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
           
            % tree overlapping group Lasso projection
            ind_work(1:2,:)=ind(1:2,:);
            ind_work(3,:)=ind(3,:) * (lambda / L);
            
            if (GFlag==0)
                x=altra_mt(v, n, k, ind_work, size(ind_work,2));
            else
                x=general_altra_mt(v, n, k, G, ind_work, size(ind_work,2));
            end
                        
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
            xRow=x(i,:);
            
            if (GFlag==0)
                tree_norm=treeNorm(xRow, k, ind, size(ind,2));
            else
                tree_norm=general_treeNorm(xRow, k, G, ind, size(ind,2));
            end
            
            funVal(iterStep)=funVal(iterStep)+ lambda*tree_norm;
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
else        
    error('\n The function does not support opts.mFlag neq 0 & opts.lFlag neq 0!');
end