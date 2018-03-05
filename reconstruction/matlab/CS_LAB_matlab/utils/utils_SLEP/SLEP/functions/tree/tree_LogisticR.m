function [x, c, funVal, ValueL]=tree_LogisticR(A, y, z, opts)
%
%%
% Function tree_LogisticR
%      Logistic Loss with the 
%           tree structured group Lasso Regularization
%
%% Problem
%
%  min  f(x,c) = - sum_i weight_i * log (p_i) + z * sum_j w_j ||x_{G_j}||
%
%  a_i denotes a training sample,
%      and a_i' corresponds to the i-th row of the data matrix A
%
%  y_i (either 1 or -1) is the response
%     
%  p_i= 1/ (1+ exp(-y_i (x' * a_i + c) ) ) denotes the probability
%
%  G_j's are nodes with tree structure
%
%  The tree overlapping group information is contained in
%  opts.ind, which is a 3 x nodes matrix, where nodes denotes the number of
%  nodes of the tree.
%  opts.ind(1,:) contains the starting index
%  opts.ind(2,:) contains the ending index
%  opts.ind(3,:) contains the corresponding weight (w_j)
%
%  Note: 
%  1) If each element of x is a leaf node of the tree and the weight for
%  this leaf node are the same, we provide an alternative "efficient" input
%  for this kind of node, by creating a "super node" with 
%  opts.ind(1,1)=-1; opts.ind(2,1)=-1; and opts.ind(3,1)=the common weight.
%
%  2) If the features are well ordered in that, the features of the left
%  tree is always less than those of the right tree, opts.ind(1,:) and
%  opts.ind(2,:) contain the "real" starting and ending indices. That is to
%  say, x( opts.ind(1,j):opts.ind(2,j) ) denotes x_{G_j}. In this case,
%  the entries in opts.ind(1:2,:) are within 1 and n.
%
%
%  If the features are not well ordered, please use the input opts.G for
%  specifying the index so that  
%   x( opts.G ( opts.ind(1,j):opts.ind(2,j) ) ) denotes x_{G_j}.
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
%  y -        Response vector (of size mx1)
%  z -        Regularization parameter (z >=0)
%  opts-      optional inputs (default value: opts=[])
%             !!For tr_LogisticR, we require that opts.ind is specified.!!
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

if (length(y) ~=m)
    error('\n Check the length of y!\n');
end

if (z<0)
    error('\n z should be nonnegative!\n');
end

opts=sll_opts(opts); % run sll_opts to set default values (flags)

%% Detailed initialization
%% Normalization

% Please refer to the function 'sll_opts'
%                 for the definitions of mu, nu and nFlag
%
% If .nFlag =1, the input matrix A is normalized to
%                     A= ( A- repmat(mu, m,1) ) * diag(nu)^{-1}
%
% If .nFlag =2, the input matrix A is normalized to
%                     A= diag(nu)^{-1} * ( A- repmat(mu, m,1) )
%
% Such normalization is done implicitly.
%     This implicit normalization is suggested for the sparse matrix,
%                                    but not for the dense matrix.
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

% Initialize ind 
if (~isfield(opts,'ind'))
    error('\n In tree_LeastR, the field .ind should be specified');
else
    ind=opts.ind;
   
    if (size(ind,1)~=3)
        error('\n Check opts.ind');
    end
end

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


GFlag=0;
% if GFlag=1, we shall apply general_altra
if (isfield(opts,'G'))
    GFlag=1;
    
    G=opts.G;
    if (max(G) >n || max(G) <1)
        error('\n The input G is incorrect. It should be within %d and %d',1,n);
    end
end
    

%% Starting point initialization

p_flag=(y==1);                  % the indices of the postive samples
m1=sum(weight(p_flag));         % the total weight for the positive samples
m2=1-m1;                        % the total weight for the positive samples

% process the regularization parameter
if (opts.rFlag==0)
    lambda=z;
else % z here is the scaling factor lying in [0,1]
%     if (z<0 || z>1)
%         error('\n opts.rFlag=1, and z should be in [0,1]');
%     end

    computedFlag=0;
    if (isfield(opts,'lambdaMax'))
        if (opts.lambdaMax~=-1)
            lambda=z*opts.lambdaMax;
            computedFlag=1;
        end
    end
    
    if (~computedFlag)
        
        % we compute ATb for computing lambda_max, when the input z is a ratio        
        b(p_flag,1)=m2;  b(~p_flag,1)=-m1;
        b=b.*weight;
        
        % compute AT b
        if (opts.nFlag==0)
            ATb =A'*b;
        elseif (opts.nFlag==1)
            ATb= A'*b - sum(b) * mu';  ATb=ATb./nu;
        else
            invNu=b./nu;               ATb=A'*invNu-sum(invNu)*mu';
        end
        
        % compute lambda_max
        if (GFlag==0)
            lambda_max=findLambdaMax(ATb, n, ind, size(ind,2));
        else
            lambda_max=general_findLambdaMax(ATb, n, G, ind, size(ind,2));
        end
        
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
        
        % we compute ATb for computing lambda_max, when the input z is a ratio        
        b(p_flag,1)=m2;  b(~p_flag,1)=-m1;
        b=b.*weight;
        
        % compute AT b
        if (opts.nFlag==0)
            ATb =A'*b;
        elseif (opts.nFlag==1)
            ATb= A'*b - sum(b) * mu';  ATb=ATb./nu;
        else
            invNu=b./nu;               ATb=A'*invNu-sum(invNu)*mu';
        end
        
        % compute lambda_max
        if (GFlag==0)
            lambda_max=findLambdaMax(ATb, n, ind, size(ind,2));
        else
            lambda_max=general_findLambdaMax(ATb, n, G, ind, size(ind,2));
        end
        
        x=lambda_max;
        funVal=lambda_max;
        ValueL=lambda_max;
        
        return;
    end
end

% initialize a starting point
if opts.init==2
    x=zeros(n,1); c=log(m1/m2);
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
        c=log(m1/m2);
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

if (opts.mFlag==0 && opts.lFlag==0)
    
    bFlag=0; % this flag tests whether the gradient step only changes a little
    
    L=1/m; % the intial guess of the Lipschitz continuous gradient
    
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
        fun_s= weight' * ( log( exp(-bb) +  exp(aa-bb) ) + bb );
        
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
        
        % copy x and Ax to xp and Axp
        xp=x;    Axp=Ax;
        cp=c;    
        
        while (1)
            % let s walk in a step in the antigradient of s to get v
            % and then do the L1/Lq-norm regularized projection
            v=s-g/L;  c= sc- gc/L;
            
            % tree overlapping group Lasso projection
            ind_work(1:2,:)=ind(1:2,:);
            ind_work(3,:)=ind(3,:) * (lambda / L);
            
            if (GFlag==0)
                x=altra(v, n, ind_work, size(ind_work,2));
            else
                x=general_altra(v, n, G, ind_work, size(ind_work,2));
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
            
            % aa= - diag(y) * (A * x + c)
            aa=- y.*(Ax+ c);
            
            % fun_x is the logistic loss at the new approximate solution
            bb=max(aa,0);
            fun_x= weight'* ( log( exp(-bb) +  exp(aa-bb) ) + bb );
            
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
                
        % compute the regularization part
        if (GFlag==0)
            tree_norm=treeNorm(x, n, ind, size(ind,2));
        else
            tree_norm=general_treeNorm(x, n, G, ind, size(ind,2));
        end
        
        % function value = loss + regularizatioin
        funVal(iterStep)=fun_x + lambda * tree_norm;    
        
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
else
    error('\n The function does not support opts.mFlag neq 0 & opts.lFlag neq 0!');
end
