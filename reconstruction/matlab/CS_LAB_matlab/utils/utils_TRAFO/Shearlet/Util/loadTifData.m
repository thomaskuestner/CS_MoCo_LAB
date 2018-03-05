function X =  loadData(fname)
info = imfinfo(fname);
num_images = numel(info);
W=info.Width;
H=info.Height;
X=zeros(H,W,num_images,'uint8');
for k = 1:num_images
    A = imread(fname, k, 'Info', info);
    X(:,:,k)=A;
end
