function y = sinc(x)

y = pi * x;

sel = abs(y) > 1e-15;

y(sel) = sin(y(sel)) ./ y(sel);
y(~sel) = 1;