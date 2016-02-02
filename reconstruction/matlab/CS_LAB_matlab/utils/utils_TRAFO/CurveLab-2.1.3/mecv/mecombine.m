function f = mecombine(r,E)
  A = size(r,1);
  N = A/2-E;
  
  f = zeros(N,size(r,2));
  
  ridx = -N-E:-N-1;  fidx = N-E:N-1;
  f(mod(fidx,N)+1,:) = f(mod(fidx,N)+1,:) - r(mod(ridx,A)+1,:);
  ridx = -N+1:-1;  fidx = N-1:-1:1;
  f(mod(fidx,N)+1,:) = f(mod(fidx,N)+1,:) + r(mod(ridx,A)+1,:);
  ridx = 0;  fidx = 0;
  f(mod(fidx,N)+1,:) = f(mod(fidx,N)+1,:) + r(mod(ridx,A)+1,:)*sqrt(2);
  ridx = 1:N-1;  fidx = 1:N-1;
  f(mod(fidx,N)+1,:) = f(mod(fidx,N)+1,:) + r(mod(ridx,A)+1,:);
  ridx = N+1:N+E-1;  fidx = N-1:-1:N-E+1;
  f(mod(fidx,N)+1,:) = f(mod(fidx,N)+1,:) - r(mod(ridx,A)+1,:);
  
  
  
  