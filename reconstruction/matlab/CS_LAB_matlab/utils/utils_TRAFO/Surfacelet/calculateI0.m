function y = calculateI0(x)

x = x / 2.0;

y = ones(size(x));
z = ones(size(x));
k = 1;

t = ones(size(x));

EPSILON = 1e-10;

while any(t > EPSILON * y)
    z = z .* x / k;
    k = k + 1;
	t = z .* z;
	y = y + t;
end
