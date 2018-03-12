function C = fdct_usfft_r2c(C)
% fdct_usfft_r2c - transform real curvelet coefficients to complex coefficients
  nbs = length(C);
  for s=2:nbs
    nw = length(C{s});
    for w=1:nw/2
      A = C{s}{w};    B = C{s}{w+nw/2};
      C{s}{w} = 1/sqrt(2) * (A+i*B);    C{s}{w+nw/2} = 1/sqrt(2) * (A-i*B);
    end
  end
