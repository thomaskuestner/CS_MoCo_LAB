function r = mescatter(f,E)
  N = size(f,1);
  A = 2*(N+E);
  
  r = zeros(A,size(f,2));
  
  ridx = -N-E:-N-1;  fidx = N-E:N-1;
  r(mod(ridx,A)+1,:) = -f(mod(fidx,N)+1,:);
  ridx = -N+1:-1;  fidx = N-1:-1:1;
  r(mod(ridx,A)+1,:) = f(mod(fidx,N)+1,:);
  ridx = 0;  fidx = 0;
  r(mod(ridx,A)+1,:) = f(mod(fidx,N)+1,:) * sqrt(2);
  ridx = 1:N-1;  fidx = 1:N-1;
  r(mod(ridx,A)+1,:) = f(mod(fidx,N)+1,:);
  ridx = N+1:N+E-1;  fidx = N-1:-1:N-E+1;
  r(mod(ridx,A)+1,:) = -f(mod(fidx,N)+1,:);
  ridx = 0;

  
  
  