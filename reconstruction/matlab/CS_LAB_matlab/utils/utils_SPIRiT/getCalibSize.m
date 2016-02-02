function [calibSize, densComp] = getCalibSize(mask)

%  [calibSize, densComp] = getCalibSize(mask)
%
%   Givem a sampling mask, the function attempts to estimate the largest region in 
%   the center of k-space that is fully sampled in order to use for calibraion.
%   it also tries to estimate the outer acceleration factor ( this will work only
%   in uniform acceleration, not variable density)
%   
% Input:
%		mask - 2D binary array
%
% Output:	
%		calibSize - [1x2] size of the calibration area
%		densComp - density compensation function.
%
%
% (c) Michael Lustig 2009


sx = 2;
sy = 2;

xflag = 0;
yflag = 0;
while 1
	
	if xflag == 0 
		tmp = crop(mask,[sx+1,sy]);
		if sum(tmp(:)) == length(tmp(:))
			sx = sx + 1;
		else
			xflag = 1;
		end
	end

	if yflag == 0 
		tmp = crop(mask,[sx,sy+1]);
		if sum(tmp(:)) == length(tmp(:))
			sy = sy + 1;
		else
			yflag = 1;
		end
	end

	if sx == size(mask,1)
		xflag = 1;
	end
	if sy == size(mask,2)
		yflag = 1;
	end

	if (xflag == 1) & (yflag == 1)
		break
	end


end

calibSize = [sx,sy];

[x,y] = meshgrid(linspace(-1,1,size(mask,2)), linspace(-1,1,size(mask,1)));
r = sqrt(x.^2+ y.^2);
circMask = r<=1;
calibMask =  zpad(ones(calibSize),[size(mask)]);

circMask = circMask - calibMask;
R = sum(circMask(:).*mask(:))/sum(circMask(:));

densComp = 1./R*(1-calibMask) + calibMask;




