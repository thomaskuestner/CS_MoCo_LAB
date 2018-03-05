function [ x_extend ] = img_extend( x, extension_y, extension_x, extension_z )
% mirrored extension to prevent any edge distortions

% M. Fischer, April 2016

%%
x_extend = zeros(size(x,1)+extension_y,size(x,2)+extension_x,size(x,3)+extension_z);
x_extend(1:end-extension_y,1:end-extension_x,1:end-extension_z) = x;
for j = 1:(extension_z-1)
    x_extend(:,:,end-extension_z+j) = x_extend(:,:,end-extension_z);
end;
x_extend(end-extension_y+1:end,1:end-extension_x,:) = flipud(x(end-extension_y+1:end,:,:));
x_extend(1:end-extension_y,end-extension_x+1:end,:) = fliplr(x(:,end-extension_x+1:end,:));
x_extend(end-extension_y+1:end,end-extension_x+1:end,:) = rot90(x(end-extension_y+1:end,end-extension_x+1:end,:),2);

end

