function X= readeAviFile(name)
xyloObj = mmreader(name);

nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;

% Preallocate movie structure.
X= zeros(vidHeight,vidWidth,nFrames);

% Read one frame at a time.
for k = 1 : nFrames
    X(:,:,k)= rgb2gray( read(xyloObj, k));
end

% Size a figure based on the video's width and height.
