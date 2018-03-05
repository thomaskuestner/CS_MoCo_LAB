function [U, out] = TVNLR(A,b,p,q,opts)

opts.TVnorm = 1;
opts.nonneg = true;
global D Dt
[D,Dt] = defDDt;
theta = opts.theta; % Zhangjian
% problem dimension
n = p*q;


% unify implementation of A
if ~isa(A,'function_handle')
    A = @(u,mode) f_handleA(A,u,mode,opts);
end

% get or check opts
opts = TVNLR_opts(opts);

% mark important constants
mu = opts.mu;
beta = opts.beta;
tol_inn = opts.tol_inn;
tol_out = opts.tol;
gam = opts.gam;

out.PSNR = zeros(opts.maxit,1);

% check if A*A'=I
tmp = rand(length(b),1);
if norm(A(A(tmp,2),1)-tmp,1)/norm(tmp,1) < 1e-3
    opts.scale_A = false;
end
clear tmp;

% check scaling A
if opts.scale_A
    [mu,A,b] = ScaleA(n,mu,A,b,opts.consist_mu);
end

% check scaling b
if opts.scale_b
    [mu,b,scl] = Scaleb(mu,b,opts.consist_mu);
end

% calculate A'*b
Atb = A(b,2);

% initialize U, beta
muf = mu;
betaf = beta;     % final beta
[U,mu,beta] = TVNLR_init(p,q,Atb,scl,opts);    % U: p*q
if mu > muf
    mu = muf;
end
if beta > betaf
    beta = betaf;
end
muDbeta = mu/beta;

rcdU = U;
X = U;% Zhangjian

% initialize multiplers
sigmax = zeros(p,q);                       % sigmax, sigmay: p*q
sigmay = zeros(p,q);
delta = zeros(length(b),1);                % delta: m

% initialize D^T sigma + A^T delta
DtsAtd = zeros(p*q,1);
gamma = zeros(p*q,1);

% initialize out.n2re
if isfield(opts,'Ut')
    Ut = opts.Ut*scl;        %true U, just for computing the error
    nrmUt = norm(Ut,'fro');
else
    Ut = [];
end
if ~isempty(Ut)
    out.n2re = norm(U - Ut,'fro')/nrmUt;
end

% prepare for iterations
out.mus = mu; out.betas = beta;
out.res = []; out.itrs = []; out.f = []; out.obj = []; out.reer = [];
out.lam1 = []; out.lam2 = []; out.lam3 = []; out.lam4 = []; out.lam5 = [];
out.itr = Inf;
out.tau = []; out.alpha = []; out.C = []; gp = [];
out.cnt = [];

[Ux,Uy] = D(U);                   % Ux, Uy: p*q
if opts.TVnorm == 1
    Wx = max(abs(Ux) - 1/beta, 0).*sign(Ux);
    Wy = max(abs(Uy) - 1/beta, 0).*sign(Uy);
    lam1 = sum(sum(abs(Wx) + abs(Wy)));
else
    V = sqrt(Ux.*conj(Ux) + Uy.*conj(Uy));        % V: p*q
    V(V==0) = 1;
    S = max(V - 1/beta, 0)./V;        % S: p*q
    Wx = S.*Ux;                       % Wx, Wy: p*q
    Wy = S.*Uy;
    lam1 = sum(sum(sqrt(Wx.*conj(Wx) + Wy.*conj(Wy))));
end

[lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Wx,Wy,lam1,beta,mu,A,b,...
    Atb,sigmax,sigmay,delta);

g3 = (U(:)-X(:)); %Zhangjian
%lam, f: constant      g2: pq        Au: m         g: pq

% compute gradient
d = g2 + muDbeta*g - DtsAtd + theta*g3;

count = 1;
Q = 1; C = f;                     % Q, C: costant
out.f = [out.f; f]; out.C = [out.C; C];
out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5];

for ii = 1:opts.maxit
    if opts.disp
        fprintf('outer iter = %d, total iter = %d, normU = %4.2e; \n',count,ii,norm(U,'fro'));
    end
    
    % compute tau first
    
    if  0%~isempty(gp)
        dg = g - gp;                        % dg: pq
        dg2 = g2 - g2p;                     % dg2: pq
        ss = uup'*uup;                      % ss: constant
        sy = uup'*(dg2 + muDbeta*dg);       % sy: constant
        % sy = uup'*((dg2 + g2) + muDbeta*(dg + g));
        % compute BB step length
        sy = abs(sy); % Zhangjian
        tau = abs(ss/max(sy,eps));               % tau: constant
        
        fst_itr = false;
    else
        % do Steepest Descent at the 1st ieration
        %d = g2 + muDbeta*g - DtsAtd;         % d: pq
        [dx,dy] = D(reshape(d,p,q));                    %dx, dy: p*q
        dDd = norm(dx,'fro')^2 + norm(dy,'fro')^2;      % dDd: cosntant
        Ad = A(d,1);                        %Ad: m
        % compute Steepest Descent step length
        tau = abs((d'*d)/(dDd + muDbeta*Ad'*Ad));
        
        % mark the first iteration
        fst_itr = true;
    end
    
    fst_itr = false;  %Zhangjian
    %tau = 0.16;fst_itr = false;  %Zhangjian
    
    % keep the previous values
    Up = U; gp = g; g2p = g2; Aup = Au; Uxp = Ux; Uyp = Uy;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ONE-STEP GRADIENT DESCENT %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    taud = tau*d;
    U = U(:) - taud;
    % projected gradient method for nonnegtivity
    if opts.nonneg
        U = max(real(U),0);
    end
    U = reshape(U,p,q);                    % U: p*q (still)
    %figure(100);imshow(U,[]);
    [Ux,Uy] = D(U);                        % Ux, Uy: p*q
    
    [lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Wx,Wy,lam1,beta,mu,A,b,...
        Atb,sigmax,sigmay,delta);
    g3 = (U(:)-X(:)); %Zhangjian
    % Nonmonotone Line Search
    alpha = 1;
    du = U - Up;                          % du: p*q
    const = opts.c*beta*(d'*taud);
    
    % Unew = Up + alpha*(U - Up)
    cnt = 0; flag = true;
    
    % if back tracking is succeceful, then recompute
    if alpha ~= 0
        Uxbar = Ux - sigmax/beta;
        Uybar = Uy - sigmay/beta;
        if opts.TVnorm == 1
            % ONE-DIMENSIONAL SHRINKAGE STEP
            Wx = max(abs(Uxbar) - 1/beta, 0).*sign(Uxbar);
            Wy = max(abs(Uybar) - 1/beta, 0).*sign(Uybar);
        else
            % TWO-DIMENSIONAL SHRINKAGE STEP
            V = sqrt(Uxbar.*conj(Uxbar) + Uybar.*conj(Uybar));
            V(V==0) = 1;
            S = max(V - 1/beta, 0)./V;
            Wx = S.*Uxbar;
            Wy = S.*Uybar;
        end
        
        % update parameters related to Wx, Wy
        [lam1,lam2,lam4,f,g2] = update_W(beta,...
            Wx,Wy,Ux,Uy,sigmax,sigmay,lam1,lam2,lam4,f,opts.TVnorm);
    end
    
    % update reference value
    Qp = Q; Q = gam*Qp + 1; C = (gam*Qp*C + f)/Q;
    uup = U - Up; uup = uup(:);           % uup: pq
    nrmuup = norm(uup,'fro');                   % nrmuup: constant
    
    out.res = [out.res; nrmuup];
    out.f = [out.f; f]; out.C = [out.C; C]; out.cnt = [out.cnt;cnt];
    out.lam1 = [out.lam1; lam1]; out.lam2 = [out.lam2; lam2]; out.lam3 = [out.lam3; lam3];
    out.lam4 = [out.lam4; lam4]; out.lam5 = [out.lam5; lam5];
    out.tau = [out.tau; tau]; out.alpha = [out.alpha; alpha];
    
    if ~isempty(Ut)
        out.n2re = [out.n2re; norm(U - Ut,'fro')/norm(Ut,'fro')];
    end
    
    % compute gradient
    d = g2 + muDbeta*g - DtsAtd + theta*g3;
    
    nrmup = norm(Up,'fro');
    RelChg = nrmuup/nrmup;
    
    NLR_flag = 1;
    if NLR_flag
        Options.kernelratio= 3;
        Options.windowratio= 6;
        Options.verbose=0;
        Options.filterstrength=0.03;
        noise_img = double(uint8(reshape((U(:)-gamma/theta)/scl,p,q)));
        J=NLMF(noise_img/255,Options);
        X = double(uint8(J*255));
        out.PSNR(ii) = csnr(X,opts.Org,0,0);        
        X = X*scl;
    else
        X = U;
    end
    
    
    if  RelChg < tol_inn && ~fst_itr
        count = count + 1;
        RelChgOut = norm(U-rcdU,'fro')/nrmup;
        out.reer = [out.reer; RelChgOut];
        rcdU = U;
        out.obj = [out.obj; f + lam4 + lam5];
        if isempty(out.itrs)
            out.itrs = ii;
        else
            out.itrs = [out.itrs; ii - sum(out.itrs)];
        end
        
        % stop if already reached final multipliers
        if RelChgOut < tol_out || count > opts.maxcnt
            if opts.isreal
                U = real(U);
            end
            if exist('scl','var')
                U = U/scl;
            end
            out.itr = ii;
            fprintf('Number of total iterations is %d. \n',out.itr);
            return
        end
        
        % update multipliers
        [sigmax,sigmay,delta,lam4,lam5,f] = update_mlp(beta,mu, ...
            Wx,Wy,Ux,Uy,Au,b,sigmax,sigmay,delta,lam4,lam5,f);
        gamma = gamma - theta*(U(:)-X(:));% Zhangjian
        % update penality parameters for continuation scheme
        beta0 = beta;
        beta = beta*opts.rate_ctn;
        mu = mu*opts.rate_ctn;
        if beta > betaf
            beta = betaf;
        end
        if mu > muf
            mu = muf;
        end
        muDbeta = mu/beta;
        out.mus = [out.mus; mu]; out.betas = [out.betas; beta];
        
        % update function value, gradient, and relavent constant
        f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;
        DtsAtd = -(beta0/beta)*d;     % DtsAtd should be divded by new beta instead of the old one for consistency!
        
        d = g2 + muDbeta*g+ theta*g3 - DtsAtd;
        %d = 2*d;  //ZhangJian
        %initialize the constants
        gp = [];
        gam = opts.gam; Q = 1; C = f;
    end
    
end

if opts.isreal
    U = real(U);
end
if exist('scl','var')
    fprintf('Attain the maximum of iterations %d. \n',opts.maxit);
    U = U/scl;
end




function [lam2,lam3,lam4,lam5,f,g2,Au,g] = get_g(U,Ux,Uy,Wx,Wy,lam1,...
    beta,mu,A,b,Atb,sigmax,sigmay,delta)
global Dt

% A*u
Au = A(U(:),1);

% g
g = A(Au,2) - Atb;



% lam2
Vx = Ux - Wx;
Vy = Uy - Wy;
lam2 = sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy)));


% g2 = D'(Du-w)
g2 = Dt(Vx,Vy);

% lam3
Aub = Au-b;
lam3 = norm(Aub,'fro')^2;

%lam4
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));

%lam5
lam5 = delta'*Aub;

% f
f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;



function [U,lam2,lam3,lam4,lam5,f,Ux,Uy,Au,g,g2] = update_g(p,q,lam1,...
    alpha,beta,mu,Up,du,gp,dg,g2p,dg2,Aup,dAu,Wx,Wy,Uxp,dUx,Uyp,dUy,b,...
    sigmax,sigmay,delta)

g = gp + alpha*dg;
g2 = g2p + alpha*dg2;
U = Up + alpha*reshape(du,p,q);
Au = Aup + alpha*dAu;
Ux = Uxp + alpha*dUx;
Uy = Uyp + alpha*dUy;

Vx = Ux - Wx;
Vy = Uy - Wy;
lam2 = sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy)));
Aub = Au-b;
lam3 = norm(Aub,'fro')^2;
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));
lam5 = delta'*Aub;
f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4 - lam5;



function [lam1,lam2,lam4,f,g2] = update_W(beta,...
    Wx,Wy,Ux,Uy,sigmax,sigmay,lam1,lam2,lam4,f,option)
global Dt

% update parameters because Wx, Wy were updated
tmpf = f -lam1 - beta/2*lam2 + lam4;
if option == 1
    lam1 = sum(sum(abs(Wx) + abs(Wy)));
else
    lam1 = sum(sum(sqrt(Wx.^2 + Wy.^2)));
end
Vx = Ux - Wx;
Vy = Uy - Wy;
g2 = Dt(Vx,Vy);
lam2 = sum(sum(Vx.*conj(Vx) + Vy.*conj(Vy)));
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));
f = tmpf +lam1 + beta/2*lam2 - lam4;



function [sigmax,sigmay,delta,lam4,lam5,f] = update_mlp(beta,mu, ...
    Wx,Wy,Ux,Uy,Au,b,sigmax,sigmay,delta,lam4,lam5,f)

Vx = Ux - Wx;
Vy = Uy - Wy;
sigmax = sigmax - beta*Vx;
sigmay = sigmay - beta*Vy;
Aub = Au-b;
delta = delta - mu*Aub;

tmpf = f + lam4 + lam5;
lam4 = sum(sum(conj(sigmax).*Vx + conj(sigmay).*Vy));
lam5 = delta'*Aub;
f = tmpf - lam4 - lam5;




function [U,mu,beta] = TVNLR_init(p,q,Atb,scl,opts)

% initialize mu beta
if isfield(opts,'mu0')
    mu = opts.mu0;
else
    error('Initial mu is not provided.');
end
if isfield(opts,'beta0')
    beta = opts.beta0;
else
    error('Initial beta is not provided.');
end

% initialize U
[mm,nn] = size(opts.init);
if max(mm,nn) == 1
    switch opts.init
        case 0, U = zeros(p,q);
        case 1, U = reshape(Atb,p,q);
    end
else
    U = opts.init*scl;
    if mm ~= p || nn ~= q
        fprintf('Input initial guess has incompatible size! Switch to the default initial guess. \n');
        U = reshape(Atb,p,q);
    end
end

function opts = TVNLR_opts(opts)
%
% Set default options.
% Written by: Chengbo Li
%
if isfield(opts,'mu')
    if ~isscalar(opts.mu) || opts.mu <0
        error('opts.mu must be positive.');
    else if opts.mu > 2^13 || opts.mu < 2^4
            fprintf('From my experience, maybe users should choose opts.mu between 2^4 and 2^13 as a priority.\n');
        end
    end
else
    opts.mu = 2^8;
end
% mu is mainly decided by noise level. Set mu big when b is noise-free
% whereas set mu small when b is very noisy.


if isfield(opts,'beta')
    if ~isscalar(opts.beta) || opts.beta <0
        error('opts.beta must be positive.');
    else if opts.beta > 2^13 || opts.beta < 2^4
            fprintf('From my experience, maybe users should choose opts.beta  between 2^4 and 2^13 as a priority.\n');
        end
    end
else
    opts.beta = 2^5;
end


% outer loop tolerence
if isfield(opts,'tol')
    if ~isscalar(opts.tol) || opts.tol <= 0
        error('opts.tol should be a positive small number.');
    end
else
    opts.tol = 1.e-6;
end;


% inner loop tolerence
if isfield(opts,'tol_inn')
    if ~isscalar(opts.tol_inn) || opts.tol_inn <= 0
        error('opts.tol_inn should be a positive small number.');
    end
else
    opts.tol_inn = 1.e-3;
end;


if isfield(opts,'maxcnt')
    if ~isscalar(opts.maxcnt) || opts.maxcnt <= 0
        error('opts.maxcnt should be a positive integer.');
    end
else
    opts.maxcnt = 10;
end


if isfield(opts,'maxit')
    if ~isscalar(opts.maxit) || opts.maxit <= 0
        error('opts.maxit should be a positive integer.');
    end
else
    opts.maxit = 1025;
end


if isfield(opts,'init')
    if length(opts.init) ~= 1
        fprintf('User has supplied opts.init as initial guess matrix......\n');
    elseif ~isinInterval(opts.init,0,1,true) || opts.init ~= floor(opts.init)
        error('opts.init should be either 0/1 or an initial guess matrix.');
    end
else
    opts.init = 1;
end


if isfield(opts,'disp')
    if ~islogical(opts.disp)
        error('opts.disp should be true or false.');
    end
else
    opts.disp = false;
end


if isfield(opts,'scale_A')
    if ~islogical(opts.scale_A)
        error('opts.scale_A should be true or false.');
    end
else
    opts.scale_A = true;
end


if isfield(opts,'scale_b')
    if ~islogical(opts.scale_b)
        error('opts.scale_b should be true or false.');
    end
else
    opts.scale_b = true;
end


if isfield(opts,'consist_mu')
    if ~islogical(opts.consist_mu)
        error('opts.consist_mu should be true or false.');
    end
else
    opts.consist_mu = false;
end
% consist_mu decides if mu should be accordingly scaled while scaling A and
% b. Strongly recommend setting as 'false' if one try to recover a signal
% or image instead of solving an exact minimization problem.


if isfield(opts,'mu0')
    if ~isscalar(opts.mu0) || opts.mu0 <= 0
        error('opts.mu0 is should be a positive number which is no bigger than beta.');
    end
else
    opts.mu0 = opts.mu;
end
% initial mu


if isfield(opts,'beta0')
    if ~isscalar(opts.beta0) || opts.beta0 <= 0
        error('opts.beta0 is should be a positive number which is no bigger than beta.');
    end
else
    opts.beta0 = opts.beta;
end
% initial beta


if isfield(opts,'rate_ctn')
    if ~isscalar(opts.rate_ctn) || opts.rate_ctn <= 1
        error('opts.rate_ctn is either not a scalar or no bigger than one.');
    end
else
    opts.rate_ctn = 2;
end
% continuation parameter for both mu and beta


if isfield(opts,'c')
    if ~isscalar(opts.c) || opts.c <= 0 || opts.c > 1
        error('opts.c should be a scalar between 0 and 1.');
    end
else
    opts.c = 1.e-5;
end


if isfield(opts,'gamma')
    if ~isscalar(opts.gamma) || opts.gamma <= 0 || opts.gamma > 1
        error('opts.gamma should be a scalar between 0 and 1.');
    end
else
    opts.gamma = .6;
end


if isfield(opts,'gam')
    if ~isscalar(opts.gam) || opts.gam <= 0 || opts.gam > 1
        error('opts.gam should be a scalar between 0 and 1.');
    end
else
    opts.gam = .9995;
end
% Control the degree of nonmonotonicity. 0 corresponds to monotone line search.
% The best convergence is obtained by using values closer to 1 when the iterates
% are far from the optimum, and using values closer to 0 when near an optimum.


if isfield(opts,'rate_gam')
    if ~isscalar(opts.rate_gam) || opts.rate_gam <= 0 || opts.rate_gam > 1
        error('opts.rate_gam should be a scalar between 0 and 1.');
    end
else
    opts.rate_gam = .9;
end
% shrinkage rate of gam


if isfield(opts,'TVnorm')
    if opts.TVnorm ~= 1 && opts.TVnorm ~= 2
        error('opts.TVnorm should be either 1(TV/L1 model) or 2(TV/L2 model).');
    end
else
    opts.TVnorm = 2;
end


if isfield(opts,'nonneg')
    if ~islogical(opts.nonneg)
        error('opts.nonneg should be true or false.');
    end
else
    opts.nonneg = false;
end


if isfield(opts,'isreal')
    if ~islogical(opts.isreal)
        error('opts.isreal should be true or false.');
    end
else
    opts.isreal = false;
end


if isfield(opts,'TVL2')
    if ~islogical(opts.TVL2)
        error('opts.TVL2 should be true or false.');
    end
else
    opts.TVL2 = false;
end
% Decide the model: TV or TV/L2. The default is TV model, which is recommended.


if isfield(opts,'tau')
    if ~isscalar(opts.tau) || opts.tau <= 0
        error('opts.tau is not positive scalar.');
    end
else
    opts.tau = 1.8;
end
function [mu,A,b] = ScaleA(n,mu,A,b,option)

% Scales mu, A and f so that the largest eigenvalue of A'*A is 1 and the
% new problem
%
% min sum_i (||wi|| + beta/2 ||Diu - wi||^2) + mu/2 ||Au - b||^2
%
% is equivalent to the old one.
%
% If option is assigned, mu will be scaled accordingly.
%

eopts.disp = 0;
eopts.tol = .05;
if ~isreal(A(rand(n,1),1))
    eopts.isreal = false;
end

fh = @(x) A(A(x,1),2);
s2 = eigs(fh,n,1,'lm',eopts);
if real(s2) > 1 + 1e-10
    if option
        mu = mu*s2;
    end
    b = b/sqrt(s2);
    A = @(x,mode) A(x,mode)/sqrt(s2);
end

return

function [mu,b,scl] = Scaleb(mu,b,option)

% Scales mu and f so that the finite difference of f is neither too small
% nor too large.
%
% If option is assigned, mu will be scaled accordingly.
%


threshold1 = .5;      % threshold is chosen by experience.
threshold2 = 1.5;
scl = 1;
b_dif = abs(max(b) - min(b));

if b_dif < threshold1
    scl = threshold1/b_dif;
    b = scl*b;
    if option
        mu = mu/scl;
    end
else if b_dif > threshold2
        scl = threshold2/b_dif;
        b = scl*b;
        if option
            mu = mu/scl;
        end
    end
end

return

