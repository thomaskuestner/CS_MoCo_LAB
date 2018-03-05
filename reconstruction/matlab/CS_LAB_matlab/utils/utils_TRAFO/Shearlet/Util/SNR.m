function s =SNR(x,xrec)
signalPower= mean(x(:).^2);
noisePower = mean((x(:) - xrec(:)).^2);
s =  signalPower/noisePower;