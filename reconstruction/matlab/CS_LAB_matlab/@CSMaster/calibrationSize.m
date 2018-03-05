function calibSize = calibrationSize(obj,mask,outlierRate,initValues)
% return the calibration size of the input k-space region
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------
    
sx = 2;
sy = 2;
sz = 2;

if(nargin < 3)
    outlierRate = [0 0 0];
end  

if(nargin == 4)
    sx = initValues(1);
    sy = initValues(2);
    sz = initValues(3);
end

if(sx >= size(mask,1))
    sx = size(mask,1); % for guaranteed correct cropping
    xflag = true;
else
    xflag = false;
end
if(sy >= size(mask,2))
    sy = size(mask,2);
    yflag = true;
else
    yflag = false;
end    
if(sz >= size(mask,3))
    sz = size(mask,3);
    zflag = true;
else
    zflag = false;
end

while(~(xflag && yflag && zflag))

    if(~xflag)
        mask_helper = crop(mask,[sx+1,sy,sz]);
        if(all(mask_helper))
            sx = sx + 1;
        else
            if(outlierRate(1) > 0)
                if(nnz(prod(prod(mask_helper,2),3) == 0) < outlierRate(1) * size(mask_helper,1))
                    sx = sx + 1;
                else
                    xflag = true;
                end
            else
                xflag = true;
            end
        end
    end

    if(~yflag)
        mask_helper = crop(mask,[sx,sy+1,sz]);
        if(all(mask_helper))
            sy = sy + 1;
        else
            if(outlierRate(2) > 0)
                if(nnz(prod(prod(mask_helper,1),3) == 0) < outlierRate(2) * size(mask_helper,2))
                    sy = sy + 1;
                else
                    yflag = true;
                end
            else
                yflag = true;
            end
        end
    end
    
    if(~zflag)
        mask_helper = crop(mask,[sx,sy,sz+1]);
        if(all(mask_helper))
            sz = sz + 1;
        else
            if(outlierRate(3) > 0)
                if(nnz(prod(prod(mask_helper,1),2) == 0) < outlierRate(3) * size(mask_helper,3))
                    sz = sz + 1;
                else
                    zflag = true;
                end
            else
                zflag = true;
            end
        end
    end

    if(sx == size(mask,1))
        xflag = true;
    end
    if(sy == size(mask,2))
        yflag = true;
    end
    if(sz == size(mask,3))
        zflag = true;
    end
end


if(size(mask,3) > 1)
    calibSize = [sx,sy,sz];
else
    calibSize = [sx,sy];
end
if(any(outlierRate < 0))
    maskSize = size(mask);
    calibSize(outlierRate < 0) = maskSize(outlierRate < 0);
end

end

