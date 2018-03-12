function w = calculateKaiserwindow(N, beta)

denominator = calculateI0(beta);

val = 2 * [0 : (N - 1) / 2] / (N - 1) - 1;
w_half = calculateI0(beta * sqrt (1 - val .* val)) / denominator;

w = [w_half w_half(end - 1 : -1 : 1)].';

