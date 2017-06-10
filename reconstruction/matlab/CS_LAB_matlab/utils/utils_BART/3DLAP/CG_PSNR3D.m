function PSNR = CG_PSNR3D(Base_image, Recon_image)

[m,n,p] = size(Recon_image);

Max_I = max(max(max(abs(Base_image))));

MSE = sum(sum(sum( (abs(Base_image-Recon_image).^2) ) ))/(n*m*p);

PSNR = 10*log10( (Max_I)^2/MSE );
return