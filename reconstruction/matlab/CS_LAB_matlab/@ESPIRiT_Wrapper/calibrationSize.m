function calibSize = calibrationSize(obj,mask)
% return the calibration size of the input k-space region
%
% (c) Thomas Kuestner 
% ---------------------------------------------------------------------

sx = 2;
sy = 2;

xflag = false;
yflag = false;
if(size(mask,3) > 1)
    zflag = false;
    sz = 2;
else
    zflag = true;
    sz = 1;
end

while(~(xflag && yflag && zflag))

    if(~xflag)
        mask_helper = crop(mask,[sx+1,sy,sz]);
        if(all(mask_helper))
            sx = sx + 1;
        else
            xflag = true;
        end
    end

    if(~yflag)
        mask_helper = crop(mask,[sx,sy+1,sz]);
        if(all(mask_helper))
            sy = sy + 1;
        else
            yflag = true;
        end
    end
    
    if(~zflag)
        mask_helper = crop(mask,[sx,sy,sz+1]);
        if(all(mask_helper))
            sz = sz + 1;
        else
            zflag = true;
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

end

