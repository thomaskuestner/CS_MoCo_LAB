function y = idct2(x)
  t = idct(x);
  t = idct(t');
  y = t';