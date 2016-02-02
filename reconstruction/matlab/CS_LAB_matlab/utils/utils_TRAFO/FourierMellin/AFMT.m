classdef AFMT < handle
% helper class for Fourier-Mellin transformation
%
% (c) Lukas Mauch, Thomas Kuestner 
% ---------------------------------------------------------------------
    
    properties
        sizeX
        sizeY
        maxSize
        maxR
        indexFW
        indexBW
        sizeLR
        sizePhi
        sigma
        extrapVal
        fwX
        fwY
        bwR
        bwPhi
        scaleFactor
    end
    
    methods
        function obj = AFMT(sizeX, sizeY, sizeLR, sizePhi, sigma, extrapVal)
            obj.sizeX       = sizeX;
            obj.sizeY       = sizeY;
            obj.sizeLR      = sizeLR;
            obj.sizePhi     = sizePhi;
            obj.sigma       = sigma;
            obj.extrapVal   = extrapVal;
            obj.maxSize     = max(sizeX,sizeY);                   
            oX              = fix(sizeX/2)+1;
            oY              = fix(sizeY/2)+1;
            obj.maxR        = sqrt(max(abs(oX-1), abs(sizeX-oX))^2 + max(abs(oY-1), abs(sizeY-oY))^2);
            
            %generate all matrices for xy to log-pol and log-pol to xy
            %conversion:                     
            [gridCartX, gridCartY]  = meshgrid( -oX+1:sizeX-oX, -oY+1:sizeY-oY );
            dR                      = log(obj.maxR)/sizeLR;
            dPhi                    = 2*pi/sizePhi;
            [gridR, gridPhi]        = meshgrid( log(obj.maxR)/sizeLR:dR:log(obj.maxR), ...
                                      linspace(0,2*pi,sizePhi) );
            gridCartR               = sqrt(gridCartX.^2 + gridCartY.^2);
            gridCartPhi = atan(gridCartY./gridCartX);
            %  unwrap:
            qA = logical(logical((gridCartX >= 0)) .* logical((gridCartY<=0)));
            qB = logical(logical((gridCartX < 0)) .* logical((gridCartY<=0)));
            qC = logical(logical((gridCartX < 0)) .* logical((gridCartY>0)));
            qD = logical(logical((gridCartX >= 0)) .* logical((gridCartY>0)));
            gridCartPhi(qA) = abs(gridCartPhi(qA));
            gridCartPhi(qB) = pi-gridCartPhi(qB);
            gridCartPhi(qC) = pi+abs(gridCartPhi(qC));
            gridCartPhi(qD) = 2*pi-gridCartPhi(qD);
            gridCartPhi(isnan(gridCartPhi)) = 0;
            
            obj.fwX     = round(exp(gridR).*cos(gridPhi) + oX);
            obj.fwX(obj.fwX<1) = 1;
            obj.fwX(obj.fwX>200) = 200;
            obj.fwY     = round(exp(gridR).*sin(gridPhi) + oY);
            obj.fwY(obj.fwY<1) = 1;
            obj.fwY(obj.fwY>200) = 200;
            obj.bwR     = round(log(gridCartR)./dR);
            obj.bwR(isinf(obj.bwR)) = 1;
            obj.bwR(obj.bwR==0)     = 1;
            
            obj.bwPhi   = round(gridCartPhi./dPhi);
            obj.bwPhi(obj.bwPhi==0) = 1;

            obj.indexFW = obj.fwX.*obj.fwY;
            obj.indexBW = obj.bwR.*obj.bwPhi;
            
            obj.scaleFactor = exp(2*pi*(gridR./max(max(gridR))) * sigma);
        end
        
        function result = fAFMT(obj, img)
            imgLogPol   = img(obj.indexFW);
            figure,imshow(uint8(abs(imgLogPol))) 
            
            result      = fftshift(fft2(obj.scaleFactor .* imgLogPol));
        end
        
    end
    
end

