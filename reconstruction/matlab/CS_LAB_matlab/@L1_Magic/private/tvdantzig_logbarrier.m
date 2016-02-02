% tvdantzig_logbarrier.m
%
% Solve the total variation Dantzig program
%
% min_x TV(x)  subject to  ||A'(Ax-b)||_\infty <= epsilon
%
% Recast as the SOCP
% min sum(t) s.t.  ||D_{ij}x||_2 <= t,  i,j=1,...,n
%                  <a_{ij},Ax - b> <= epsilon  i,j=1,...,n
% and use a log barrier algorithm.
%
% Usage:  xp = tvdantzig_logbarrier(x0, A, At, b, epsilon, lbtol, mu, cgtol, cgmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a K 
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
%
% b - Kx1 vector of observations.
%
% epsilon - scalar, constraint relaxation parameter
%
% lbtol - The log barrier algorithm terminates when the duality gap <= lbtol.
%         Also, the number of log barrier iterations is completely
%         determined by lbtol.
%         Default = 1e-3.
%
% mu - Factor by which to increase the barrier constant at each iteration.
%      Default = 10.
%
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%     Default = 1e-8.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = tvdantzig_logbarrier(x0, A, At, b, epsilon, lbtol, mu, cgtol, cgmaxiter, bComplex, sz)  

largescale = isa(A,'function_handle');

if (nargin < 6), lbtol = 1e-3; end
if (nargin < 7), mu = 10; end
if (nargin < 8), cgtol = 1e-8; end
if (nargin < 9), cgmaxiter = 200; end

newtontol = lbtol;
newtonmaxiter = 50;

N = length(x0);
% n = round(sqrt(N));
if(nargin > 10 && exist('sz','var'))
    n = sz(1);
    m = sz(2);
    if(length(sz) == 3)
        o = sz(3);
        b3D = true;
    else
        o = 1;
        b3D = false;
    end
    
%     if(bComplex)
%         N = 2*n*m;
%     else
%         N = n*m;
%     end
else
    % quadratic FOV (assume 2D!)
    b3D = false;
    if(bComplex)
        n = round(sqrt(N/2));
        m = n;
    else
        n = round(sqrt(N));
        m = n;
    end
    sz = [n, m];
end

% create (sparse) differencing matrices for TV
% Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
%   reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
% Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
%   reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);
if(bComplex)    
    if(b3D)
        Dv = spdiags([repmat([reshape([-ones(n-1,m); zeros(1,m)],m*n,1) ...
          reshape([zeros(1,m); ones(n-1,m)],m*n,1)],o,1); repmat([reshape([-ones(n-1,m); zeros(1,m)],m*n,1) ...
          reshape([zeros(1,m); ones(n-1,m)],m*n,1)],o,1)], [0 1], N, N);
        Dh = spdiags([repmat([reshape([-ones(n,m-1) zeros(n,1)],m*n,1) ...
          reshape([zeros(n,1) ones(n,m-1)],m*n,1)],o,1); repmat([reshape([-ones(n,m-1) zeros(n,1)],m*n,1) ...
          reshape([zeros(n,1) ones(n,m-1)],m*n,1)],o,1)], [0 n], N, N);
        Dd = spdiags([repmat([-ones(m*n,1) ones(m*n,1)],o-1,1); [zeros(m*n,1), ones(m*n,1)]; ...
            [-ones(m*n,1) zeros(m*n,1)]; repmat([-ones(m*n,1) ones(m*n,1)],o-2,1); [zeros(m*n,1), ones(m*n,1)]], [0, n*m], N, N);      

%         Dd = spdiags([repmat([-ones(m*n,1) ones(m*n,1)],o,1); repmat([-ones(m*n,1) ones(m*n,1)],o,1)], [0, n*m], N, N);
%         Dd(m*n*(o-1)+1:m*n*o,m*n*(o-1)+1:m*n*o) = 0; % last slice
%         Dd(2*m*n*(o-1)+1:2*m*n*o,2*m*n*(o-1)+1:2*m*n*o) = 0;
%         Dd(1:m*n*o,m*n*o+1:end) = 0; % no interaction between real and imaginary part
%         Dd(m*n*o+1:end,1:m*n*o) = 0;
        
%         Dd(end+1:end+N/2,end+1:end+N/2) = Dd;
    else
        Dv = spdiags([reshape([-ones(n-1,m); zeros(1,m)],N/2,1) ...
          reshape([zeros(1,m); ones(n-1,m)],N/2,1); reshape([-ones(n-1,m); zeros(1,m)],N/2,1) ...
          reshape([zeros(1,m); ones(n-1,m)],N/2,1)], [0 1], N, N);
        Dh = spdiags([reshape([-ones(n,m-1) zeros(n,1)],N/2,1) ...
          reshape([zeros(n,1) ones(n,m-1)],N/2,1); reshape([-ones(n,m-1) zeros(n,1)],N/2,1) ...
          reshape([zeros(n,1) ones(n,m-1)],N/2,1)], [0 n], N, N);
    end
else
    if(b3D)
        Dv = spdiags(repmat([reshape([-ones(n-1,m); zeros(1,m)],m*n,1) ...
            reshape([zeros(1,m); ones(n-1,m)],m*n,1)],o,1), [0 1], N, N);
        Dh = spdiags(repmat([reshape([-ones(n,m-1) zeros(n,1)],m*n,1) ...
            reshape([zeros(n,1) ones(n,m-1)],m*n,1)],o,1), [0 n], N, N);
        Dd = spdiags([repmat([-ones(m*n,1) ones(m*n,1)],o-1,1); [zeros(m*n,1), ones(m*n,1)]], [0, n*m], N, N); 
%         Dd = spdiags(repmat([-ones(m*n,1) ones(m*n,1)],o,1), [0, n*m], N, N);
%         Dd(m*n*(o-1)+1:end,m*n*(o-1)+1:end) = 0;
%         Dd = spdiags([-ones(N-m*n,1) ones(N-m*n,1); zeros(m*n,1) ones(m*n,1)], [0, n*m], N, N);
    else
        Dv = spdiags([reshape([-ones(n-1,m); zeros(1,m)],N,1) ...
          reshape([zeros(1,m); ones(n-1,m)],N,1)], [0 1], N, N);
        Dh = spdiags([reshape([-ones(n,m-1) zeros(n,1)],N,1) ...
          reshape([zeros(n,1) ones(n,m-1)],N,1)], [0 n], N, N);
    end
end

if (largescale)
  if (norm(A(x0)-b) > epsilon)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    AAt = @(z) A(At(z));
    [w, cgres] = cgsolve(AAt, b, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = At(w);
  end
else
  if (norm(A*x0-b) > epsilon)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    opts.POSDEF = true; opts.SYM = true;
    [w, hcond] = linsolve(A*A', b, opts);
    if (hcond < 1e-14)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = A'*w;
  end  
end
x = x0;
Dhx = Dh*x;  Dvx = Dv*x;
if(b3D), Ddx = Dd*x; end;
if(bComplex)
    tmph = Dhx(1:N/2) + 1i * Dhx(N/2+1:N);
    tmpv = Dvx(1:N/2) + 1i * Dvx(N/2+1:N);
    if(b3D), tmpd = Ddx(1:N/2) + 1i * Ddx(N/2+1:N); end;
end
if(b3D)
    t = 1.05*sqrt(Dhx.^2 + Dvx.^2 + Ddx.^2) + .01*max(sqrt(Dhx.^2 + Dvx.^2 + Ddx.^2));
else
    t = 1.05*sqrt(Dhx.^2 + Dvx.^2) + .01*max(sqrt(Dhx.^2 + Dvx.^2));
end

% choose initial value of tau so that the duality gap after the first
% step will be about the origial TV
if(bComplex)
    if(b3D)
        tau = 3*(N/2)/sum(sqrt(abs(tmph).^2 + abs(tmpv).^2 + abs(tmpd).^2));
    else
        tau = 3*(N/2)/sum(sqrt(abs(tmph).^2 + abs(tmpv).^2));
    end
    lbiter = ceil((log(3*(N/2))-log(lbtol)-log(tau))/log(mu));
else
    if(b3D)
        tau = 3*N/sum(sqrt(Dhx.^2+Dvx.^2+Ddx.^2));
    else
        tau = 3*N/sum(sqrt(Dhx.^2+Dvx.^2));
    end
    lbiter = ceil((log(3*N)-log(lbtol)-log(tau))/log(mu));
end
    
% disp(sprintf('Number of log barrier iterations = %d\n', lbiter));
totaliter = 0;
dispProgress('Log barrier', 0, lbiter);
for ii = 1:lbiter
  
  [xp, tp, ntiter] = tvdantzig_newton(x, t, A, At, b, epsilon, tau, newtontol, newtonmaxiter, cgtol, cgmaxiter, bComplex, b3D, sz);
  totaliter = totaliter + ntiter;
%   tvxp = sum(sqrt((Dh*xp).^2 + (Dv*xp).^2));
  
%   disp(sprintf('\nLog barrier iter = %d, TV = %.3f, functional = %8.3f, tau = %8.3e, total newton iter = %d\n', ...
%     ii, tvxp, sum(tp), tau, totaliter));
  dispProgress('Log barrier',ii/lbiter);

  x = xp;
  t = tp;
  
  tau = mu*tau;
  
end
dispProgress('Log barrier', 'Close');
end
                   