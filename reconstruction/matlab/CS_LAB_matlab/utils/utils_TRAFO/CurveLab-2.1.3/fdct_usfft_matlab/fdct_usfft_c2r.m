function C = fdct_usfft_c2r(C)
% fdct_usfft_c2r - transform complex curvelet coefficients to real coefficients
  nbs = length(C);
  C{1}{1} = real(C{1}{1});
  for s=2:nbs
    nw = length(C{s});
    for w=1:nw/2
      A = C{s}{w};      %B = C{s}{w+nw/2};
      C{s}{w} = sqrt(2) * real(A);      C{s}{w+nw/2} = sqrt(2) * imag(A);
    end
  end
  C{nbs}{1} = real(C{nbs}{1});
    
  
