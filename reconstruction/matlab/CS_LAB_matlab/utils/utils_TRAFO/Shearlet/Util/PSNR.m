function PSNRdb = PSNR(x, y)

err = x - y;
err = err(:);
PSNRdb = 20 * log10(256/sqrt(mean(err .^2)));