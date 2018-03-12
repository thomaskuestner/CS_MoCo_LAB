% April 7 2009
% written by T. K. Pong and P. Tseng
%
% This Matlab function solves
%   min_{W} f(W) = 1/2*\|A*W-Y\|^2_F+\|W\|_*,
% where A is p by n and Y is p by m (so W is n by m), and
% \|\|_* is the nuclear norm.
% The code requires that p>=m.
% 
% This code first reduces A to have full row rank.  Then it
% applies an accelerated proximal gradient method
% to solve the primal problem.
%
% This matlab function can be called using the command
%
% [W,fmin] = mat_primal(A,Y,lambda,opt);
%
% Input:
% A,Y are as described above.
% lambda: regularization parameter
% opt.tol    = termination tolerance (e.g., 1e-3).
% opt.freq   = frequency of termination checks (e.g., n or 10).
%
% Output: 
% W  = the optimal solution
% pobj   = optimal objective function value of the primal problem.

function [W,pobj]=mat_primal(A,Y,lambda,opt)

B = Y;
clear Y;

if nargin<4
    opt = [];
end

if isfield(opt, 'tol')
    tol = opt.tol;
else
    tol = 10^-3;
end

if isfield(opt, 'freq')
    freq = opt.freq;
else
    freq = 10;
end

A = (1/sqrt(lambda))*A;
B = (1/sqrt(lambda))*B;

n=size(A,2);
m=size(B,2);
p=size(B,1);
fprintf(' m,n,p = %g,%g,%g\n', m,n,p)
maxiter=1000000;
meps=1e-10;
iter=0;

% check if A has full row rank
reduce_flag=1;
if p>=n 
  r=rank(A);		% rank is expensive, do this only if needed
  if r==n 
    reduce_flag=0;
  end
end

% If A lacks full row rank, reduce A to "upper triangular".  
if reduce_flag==1
  fprintf(' reduce A to have full row rank: \n');
  tic
  [R0,S0,E0] = qr(A');
  r=rank(S0);
  Anew=(S0(1:r,:))';
  Bnew=E0'*B;
  clear S0 E0
  t_rA=toc;
  fprintf(' done reducing A, time: %g\n', t_rA);
else
  r=n;
  t_rA=0;
end
 


% From now on, the n in the code and in the description corresponds to the
% rank of Anew, i.e., r.
tic

if reduce_flag==1 
  M=Anew'*Anew;
  M2=Anew'*Bnew;
  m3=mytrace(Bnew,Bnew);
%  C=inv(M);
%  E=Bnew'*Anew*C;
else
  M=A'*A;
  M2=A'*B;
  m3=mytrace(B,B);
%  C=inv(M);
%  E=B'*A*C;
end
C=inv(M);
E=M2'*C;
f=mytrace(M2*M2',C)-m3;

t3=toc;
fprintf(' done computing C and E, time: %g \n',t3);

tic
% Initialize W
W=E';
%W=C*M2;  	%This least-square init yields comparable performance as W=E'.
%W=zeros(r,m);	%This yields much slower convergence.
W0=W;

theta=1;
theta0=1;

options.disp=0;
L=eigs(M,1,'LM',options); 	%Lipschitz constant (alg diverges if it's halved)
%L=norm(M);		%slower than eigs(M,1)
fprintf(' L = %g , tol = %g, freq = %g \n', L, tol, freq)

while iter<=maxiter
  
  iter=iter+1;

  Y=W+theta*(1/theta0-1)*(W-W0);
  G=M*Y-M2;
  T=Y-G/L;
  [R,D,S]=svd(T,'econ'); % T=R*D*S'; This and the following line compute
                           % the minimizer to step 2 of accel. grad.
                           % algorithm.
  W0=W;
  W=R*(max(D-eye(size(D,1))/L,0))*S';
%  W=R*diag(max(diag(D)-1/L,0))*S';
%  W=R*diag(median([diag(D)+1/L,diag(D)-1/L,zeros(min(m,n),1)]'))*S';

  if norm(W-W0,'fro')<=meps
    fprintf(' termination due to negligible change in U = %g \n',norm(W-W0,'fro'));
    U=(W'-E)*M;
    [R,D,S]=svd(U,'econ');
    U=R*min(D,1)*S'; 		%Project to make U dual feasible
    pobj=(mytrace(W,M*W)-2*mytrace(M2,W)+m3)/2+sum(svd(W));
    dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
    fprintf(' iter= %g  dobj= %g  pobj= %g \n',iter,dobj,pobj);
    break
  end
  
  theta0=theta;
  theta=(sqrt(theta^4+4*theta^2)-theta^2)/2; % Update theta
  iter=iter+1;

  % Compute the gradient of the smooth part
  if (iter>0) && (mod(iter,freq)==0)
    U=(W'-E)*M;
    [R,D,S]=svd(U,'econ');
    U=R*min(D,1)*S';		%Project to make U dual feasible
    pobj=(mytrace(W,M*W)-2*mytrace(M2,W)+m3)/2+sum(svd(W));
    dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
    fprintf(' iter= %g  dobj= %g  pobj= %g \n',iter,dobj,pobj);
    if abs(pobj-dobj) < tol*(abs(dobj)+1)
      break
    end
  end
end
t1=toc;
if reduce_flag==1
  W=R0*[W;zeros(n-r,m)];
end
fprintf(' iter = %g, fmin = %g, total time = %g\n', iter, pobj, t1+t3+t_rA)