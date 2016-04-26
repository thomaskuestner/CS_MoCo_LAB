function sd_estimate = sdest3(x) 
 wc=dwt3(x,'sym4','mode','sym');
s=wc.dec{2,2,2};
sd_estimate = mad(s(:),1)/0.6745;
