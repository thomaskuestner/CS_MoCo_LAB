% revised on September 21, 2010 by Shuiwang Ji to make parameters consistent

% March 29 2009
% written by T. K. Pong and P. Tseng
%
% This Matlab function solves
%   min_{W} f(W) = 1/2*\|A*W-Y\|^2_F+\|W\|_*,
% where A is p by n and Y is p by m (so W is n by m), and
% \|\|_* is the nuclear norm.
% The code requires that p>=m.
% 
% This code first reduces A to have full row rank.  Then it
% applies a feasible descent method to solve the dual
% problem:
%   min_U 1/2*<C,U*U'> + <E',U>  s.t.   U'*U <= I, 
% where C=(A'*A)^{-1} and E=Y'*A*C (so U is m by n).
% From the solution U we obtain  W = E'+C*U'.
%
% This matlab function can be called using the command
%
% [W,fmin] = mat_dual(A,Y,lambda,opt);
%
% Input:
% A,Y are as described above.
% lambda: regularization parameter
% opt.alg    = 1  Frank-Wolfe method with line search 
%              2  gradient-projection method with large constant stepsize + line search
%              3  Accelerated grad-proj method with small constant stepsize
% opt.tol    = termination tolerance (e.g., 1e-3).
% opt.freq   = frequency of termination checks (e.g., n or 10).
%
% Output: 
% W  = the optimal solution
% pobj   = optimal objective function value of the original problem.

function [W,pobj] = mat_dual(A,Y,lambda,opt)

B = Y;
clear Y;

if nargin<4
    opt = [];
end

if isfield(opt, 'alg')
    alg = opt.alg;
else
    alg = 3;
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

% Read data parameters
n=size(A,2);
m=size(B,2);
p=size(A,1);
if size(B,1)~=p
  error('A and B must have same number of rows');
end

% termination tolerance & frequency of termination checks
%tol=input(' enter the tolerance (1e-3) ');
fprintf(' n = %g  m = %g  p= %g  tol= %g \n',n,m,p,tol);
maxiter=1000000;		% maximum number of iterations
meps=1e-10;			% machine epsilon: used for early termination check
                        % & perturbing for a nearly optimal solution 

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
  nold=n;
  n=r;
  t_rA=toc;
  fprintf(' done reducing A, time: %g\n', t_rA);
else
  t_rA=0;
end
 


% From now on, the n in the code and in the description corresponds to the
% rank of Anew, i.e., r.
tic


if reduce_flag==1 
  M=Anew'*Anew;
  M2=2*Anew'*Bnew;
  m3=mytrace(Bnew,Bnew);
  C=inv(M);
  E=Bnew'*Anew*C;
  f=mytrace(Bnew*E,Anew)-m3;
else
  M=A'*A;
  M2=2*A'*B;
  m3=mytrace(B,B);
  C=inv(M);
  E=B'*A*C;
  f=mytrace(B*E,A)-m3;
end
t0=toc;
fprintf(' done computing C and E, time: %g \n',t0);

  if alg==1
    k=0;
    stop=0;
    tic

    U=zeros(m,n);    % initialize with zero U.
    options.disp=0;  % suppress display in eigs
    
    while stop==0 && k <=maxiter
      G=U*C+E;
      [R,D,S] = svd(G,'econ');    % For m<=n, R and D are m by m, S is n by m; otherwise, S and D are n by n, R is m by n.
      hU=-R*S';
      Del=hU-U;           % Del is the search direction
      alpha=min(1,-mytrace(G,Del)/mytrace(Del*C,Del));
      
      % early termination if stepsize goes negative
      if alpha<meps
        fprintf(' termination due to negative stepsize = %g \n',alpha);
        W=E'+ C*U';
        pobj=(mytrace(W,M*W)-mytrace(M2,W)+m3)/2+sum(svd(W));
        dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
        dfeas=max(0,svds(U,1)^2-1);
        fprintf(' iter= %g  dobj= %g  dual feas= %g  pobj=  %g \n',k,dobj,dfeas,pobj);
        break
      end
      
      U=U+alpha*Del;
      k=k+1;    
      if (k>0) && (mod(k,freq)==0)
        W=E'+ C*U';
        pobj=(mytrace(W,M*W)-mytrace(M2,W)+m3)/2+sum(svd(W));
        dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
        dfeas=max(0,svds(U,1)^2-1);
        fprintf(' iter= %g  dobj= %g  dual feas= %g  pobj= %g\n',k,dobj,dfeas,pobj);
        if abs(pobj-dobj) < tol*(abs(dobj)+1) && dfeas<tol
          stop=1;
        end
      end
    end
    t1=toc;
    fprintf(' Frank-Wolfe: iter= %g  total_time= %g \n',k,t1+t0+t_rA);
  end
  
  if alg==2
    k=0;
    stop=0;
    tic
    
    % Do one Frank-Wolfe iteration from zero U.  This works better than zero U.
    [R,D,S] = svd(E,'econ');    % For m<=n, R and D are m by m, S is n by m; otherwise, S and D are n by n, R is m by n.
    Del=-R*S';
    alpha=min(1,-mytrace(E,Del)/mytrace(Del*C,Del));
    U=alpha*Del;
    options.disp=0;  % suppress display in eigs

%    L=norm(C)/8;
    L=eigs(C,1,'LM',options)/8;
    while stop==0 && k <=maxiter
      G=U*C+E;
      [R,D,S] = svd(U-G/L,'econ');    % For m<=n, R and D are m by m, S is n by m; otherwise, S and D are n by n, R is m by n.
      hU=R*min(D,eye(min(m,n)))*S';
      Del=hU-U;           % Del is the search direction
      alpha=min(1,-mytrace(G,Del)/mytrace(Del*C,Del));
      
      % early termination if stepsize goes negative
      if alpha<meps
        fprintf(' termination due to negative stepsize = %g \n',alpha);
        W=E'+ C*U';
        pobj=(mytrace(W,M*W)-mytrace(M2,W)+m3)/2+sum(svd(W));
        dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
        dfeas=max(0,svds(U,1)^2-1);
        fprintf(' iter= %g  dobj= %g  dual feas= %g  pobj= %g \n',k,dobj,dfeas,pobj);
        break
      end
      
      U=U+alpha*Del;
      k=k+1;  
      if (k>0) && (mod(k,freq)==0)
        W=E'+ C*U';
        pobj=(mytrace(W,M*W)-mytrace(M2,W)+m3)/2+sum(svd(W));
        dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
        dfeas=max(0,svds(U,1)^2-1);
        fprintf(' iter= %g  dobj= %g  dual feas= %g  pobj= %g \n',k,dobj,dfeas,pobj);
        if abs(pobj-dobj) < tol*(abs(dobj)+1) && dfeas<tol
          stop=1;
        end
      end
    end
    t1=toc;
    fprintf(' grad-proj+LS: iter= %g  total_time= %g\n',k,t1+t0+t_rA);
  end
  
  if alg==3
    k=0;
    stop=0;
    tic

    U=zeros(m,n);    % initialize with zero U.
    
    U0=U;
    theta=1;
    theta0=1;
    %  L=norm(C);
    options.disp=0;   % suppress display in eigs
    L=eigs(C,1,'LM',options)/2;			% works fine in practice
    fprintf(' L = %g\n', L);
    
    while stop==0 && k <=maxiter
      Y=U+(theta/theta0-1)*(U-U0);
      G=Y*C+E;
      [R,D,S] = svd(Y-G/L,'econ');    % For m<=n, R and D are m by m, S is n by m; otherwise, S and D are n by n, R is m by n.
      U0=U;
      U=R*min(D,eye(min(m,n)))*S';

      % early termination if U-U0 is tiny
      if norm(U-U0,'fro')<meps
        fprintf(' termination due to negligible change in U = %g \n',norm(U-U0,'fro'));
        W=E'+ C*U';
        pobj=(mytrace(W,M*W)-mytrace(M2,W)+m3)/2+sum(svd(W));
        dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
        dfeas=max(0,svds(U,1)^2-1);
        fprintf(' iter= %g  dobj= %g  dual feas= %g  pobj= %g \n',k,dobj,dfeas,pobj);
        break
      end

      theta0=theta;
      theta=2*theta/(sqrt(theta^2+4)+theta);
      k=k+1;    
      if (k>0) && (mod(k,freq)==0)
        W=E'+ C*U';
        pobj=(mytrace(W,M*W)-mytrace(M2,W)+m3)/2+sum(svd(W));
        dobj=-mytrace(U*C,U)/2 - mytrace(E,U)-f/2;
        dfeas=max(0,svds(U,1)^2-1);
        fprintf(' iter= %g  dobj= %g  dual feas= %g  pobj= %g \n',k,dobj,dfeas,pobj);
        if abs(pobj-dobj) < tol*(abs(dobj)+1) && dfeas<tol
          stop=1;
        end
      end
    end
    t1=toc;
    fprintf(' accel. grad-proj: iter= %g  total_time= %g\n',k,t1+t0+t_rA);
  end
  if reduce_flag==1
    W=R0*[W;zeros(nold-n,m)];
  end


                            