function funs = Filter_Functions
  funs.Filter_General=@Filtering;
  funs.Filter_Gauss=@Filtering_Gauss;
  funs.Filter_Laplacian=@Filtering_Laplacian;
end

function I_out = Filtering_Gauss(I_in, sigma, x)
% function performs gaussian filtering on a 3D object I_in. The gaussian
% filter is assumed to isotropic in nature (i.e. sigma_x = sigma_y =
% sigma_z).
%

F=exp(-x.*x/2/sigma^2);
F = F./sum(F(:));
I_out = imfilter(I_in, reshape(F,[length(F) 1 1]),'symmetric');
I_out = imfilter(shiftdim(I_out,1), reshape(F,[length(F) 1 1]),'symmetric');
I_out = imfilter(shiftdim(I_out,1), reshape(F,[length(F) 1 1]),'symmetric');
I_out = shiftdim(I_out,1);
end

function I_out = Filtering(I_in, F)
% function performs filtering on a 3D object I_in.
%
I_out = imfilter(I_in, reshape(F,[length(F) 1 1]),'symmetric');
I_out = imfilter(shiftdim(I_out,1), reshape(F,[length(F) 1 1]),'symmetric');
I_out = imfilter(shiftdim(I_out,1), reshape(F,[length(F) 1 1]),'symmetric');
I_out = shiftdim(I_out,1);
end

function I_out = Filtering_Laplacian(I_in)
% function performs Laplacian filtering on a 3D object I_in.
%

F = [sqrt(0.5), 1/sqrt(0.5), sqrt(0.5)];

I_out = imfilter(I_in, reshape(F,[length(F) 1 1]),'symmetric');
I_out = imfilter(shiftdim(I_out,1), reshape(F,[length(F) 1 1]),'symmetric');
I_out = imfilter(shiftdim(I_out,1), reshape(F,[length(F) 1 1]),'symmetric');
I_out = shiftdim(I_out,1);

I_out = I_out - (length(F).^3-1).*I_in;

end