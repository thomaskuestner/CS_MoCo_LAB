function x = dct2(y);
  t = dct(y);
  t = dct(t');
  x = t';
  