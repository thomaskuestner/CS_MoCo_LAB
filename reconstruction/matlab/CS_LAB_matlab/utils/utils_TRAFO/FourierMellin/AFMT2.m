function [ Y,gR,gPhi ] = AFMT2( X, sigma, mAFMT, nAFMT, extrapVal, interpMethod, algorithm, para )
%Computes discrete approximation to analytic fourier-mellin transformation.
%Parameters:
%           x                   ...matrix with input data
%           sigma               ...strictly positive number, ensuring
%                                  existence of transformation
%           mAFMT               ...size of mellin-transformation 
%
%           nAFMT               ...size of fourier-transformation
%           para                ...fast/accurate (type of backtrafo interpolation)

%algorithm selects: F-AFMT      ...fast AFMT
%                   D-AFMT      ...discrete AFMT
%
% (c) Lukas Mauch, Thomas Kuestner 
% ---------------------------------------------------------------------

switch algorithm
    %% forward-transformation:
    case 'F-AFMT' 
        [yLogPol, gR, gPhi] = convToLogPol(X, mAFMT, nAFMT, interpMethod);
       % imshow(uint8(abs(yLogPol)))
        scaleFactor = exp(2*pi*(gR./max(max(gR))) * sigma); %maximal radius normalized to 2*pi 
        Y = fftshift(fft2(scaleFactor .* yLogPol));
        
    case 'D-AFMT'
        
    %% backward-transformation:
    case 'iF-AFMT' 
        [sizeY, sizeX]  = size(X);
        
        %Compute center positions oM and oN:
        oY = fix(mAFMT/2)+1;
        oX = fix(nAFMT/2)+1;
        
        %Cartesian sample points:
        [gridCartX, gridCartY] = meshgrid( -oX+1:nAFMT-oX, -oY+1:mAFMT-oY );
        
        %Log-polar plane sample points:
        R  = sqrt(max(abs(oY-1), abs(mAFMT-oY))^2 + max(abs(oX-1), abs(nAFMT-oX))^2);        
        [gridR, gridPhi]   = meshgrid( log(R)/sizeX:log(R)/sizeX:log(R), ...
                                       linspace(0,2*pi,sizeY) );
        %0:2*pi/sizeY:2*pi/sizeY*(sizeY-1)
        
        scaleFactor = exp(-2*pi*sigma*(gridR./max(max(gridR))));
        yLogPol = scaleFactor .* ifft2(ifftshift(X));
        
        Y = convToCart(yLogPol, gridR, gridPhi, gridCartX, gridCartY, interpMethod);
      %  imshow(uint8(abs(Y)))
    case 'iD-AFMT'    
        
    otherwise 
    error('Unknown type of algorithm!')
    
end


    %% helper functions:
    function [xLogPol, gridR, gridPhi] = convToLogPol( xCart, sizeR, sizePhi, interpMethod )
        [m, n]  = size(xCart);
        
        %Compute center positions oM and oN:
        oM = fix(m/2)+1;
        oN = fix(n/2)+1;
        
        %Cartesian sample points:
        [gridX, gridY] = meshgrid( -oN+1:n-oN, fliplr(-oM+1:m-oM) );
        
        %Log-polar sample points:
        maxR = sqrt(max(abs(oN-1), abs(n-oN))^2 + max(abs(oM-1), abs(m-oM))^2);
        [gridR, gridPhi]   = meshgrid( log(maxR)/sizeR:log(maxR)/sizeR:log(maxR), ...
                                      linspace(0,2*pi,sizePhi) );

        %% Perform conversion to log-polar coordinates:

        %Assume periodicity of data in radial direction (with respect to inscribed circle):
%         minR = min([abs(oN-1); abs(n-oN); abs(oM-1); abs(m-oM)]); %determine radius of inscribed circle
%         convGridX   = mod(exp(gridR), minR) .* cos(gridPhi);
%         convGridY   = mod(exp(gridR), minR) .* sin(gridPhi);

	%Assume periodicity of data in radial direction (with respect to circumscribed circle):
        convGridX   = exp(gridR) .* cos(gridPhi);
        convGridY   = exp(gridR) .* sin(gridPhi);


        xLogPol     = interp2(gridX, gridY, xCart, convGridX, convGridY, interpMethod, extrapVal);
    end

    function xCart = convToCart( xLogPol, gridR, gridPhi, gridX, gridY, interpMethod )            
        if strcmp(para, 'fast')
        %interpolation in r-phi plane (worse but fast):
            convGridR   = sqrt(abs(gridX).^2 + abs(gridY).^2);
            convGridPhi = atan(gridY./gridX);
            qA = logical(logical((gridX >= 0)) .* logical((gridY<=0)));
            qB = logical(logical((gridX < 0)) .* logical((gridY<=0)));
            qC = logical(logical((gridX < 0)) .* logical((gridY>0)));
            qD = logical(logical((gridX >= 0)) .* logical((gridY>0)));
            convGridPhi(qA) = abs(convGridPhi(qA));
            convGridPhi(qB) = pi-convGridPhi(qB);
            convGridPhi(qC) = pi+abs(convGridPhi(qC));
            convGridPhi(qD) = 2*pi-convGridPhi(qD);
            convGridPhi(isnan(convGridPhi)) = 0;
            xCart       = interp2(exp(gridR), gridPhi, xLogPol, convGridR, convGridPhi, interpMethod, extrapVal);
        elseif strcmp(para, 'accurate')
            %Perform conversion to cartesian coordinates, using interpolation in xy-Plane (good but slow):
            convGridX   = reshape(exp(gridR) .* cos(gridPhi),[], 1);
            convGridY   = reshape(exp(gridR) .* sin(gridPhi),[], 1);

            F = TriScatteredInterp(convGridX, convGridY, reshape(xLogPol, [], 1), interpMethod);
            xCart       = F(gridX, flipud(gridY));
            xCart(isnan(xCart)) = 0;
        end
    end

end

