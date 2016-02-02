function [ Y, gR, gPhi, gTheta ] = AFMT3( X, sigma, mAFMT, nAFMT, pAFMT, extrapVal, interpMethod, algorithm )
%Computes discrete approximation to analytic fourier-mellin transformation.
%Parameters:
%           x                   ...matrix with input data
%           sigma               ...strictly positive number, ensuring
%                                  existence of transformation
%           mAFMT               ...size of mellin-transformation (check)
%
%           nAFMT               ...size of fourier-transformation (check)

%algorithm selects: F-AFMT      ...fast AFMT
%                   D-AFMT      ...discrete AFMT
%
% (c) Lukas Mauch, Thomas Kuestner 
% ---------------------------------------------------------------------

switch algorithm
    %% forward-transformation (3d):
    case 'F-AFMT' 
        [yLogPol, gR, gPhi, gTheta] = convToLogPol(X, mAFMT, nAFMT, pAFMT, interpMethod);
        scaleFactor = exp(2*pi*(gR./max(max(max(gR)))) * sigma); %maximal radius normalized to 2*pi 
        Y = fftnshift(scaleFactor .* yLogPol, 1:3);
        
    case 'D-AFMT'
        
    %% backward-transformation:
    case 'iF-AFMT' 
        [sizeY, sizeX, sizeZ]  = size(X);
        
        %Compute center positions oM and oN:
        oY = fix(mAFMT/2)+1;
        oX = fix(nAFMT/2)+1;
        oZ = fix(pAFMT/2)+1;
        
        %Cartesian sample points:
        [gridCartX, gridCartY, gridCartZ] = meshgrid( -oX+1:nAFMT-oX, -oY+1:mAFMT-oY, -oZ+1:pAFMT-oZ );
        
        %Log-polar plane sample points:
        R  = sqrt(max(abs(oY-1), abs(mAFMT-oY))^2 + max(abs(oX-1), abs(nAFMT-oX))^2 + ...
             max(abs(oZ-1), abs(pAFMT-oZ))^2);        
        [gridR, gridPhi, gridTheta]   = meshgrid( log(R)/sizeX:log(R)/sizeX:log(R), ...
                                       linspace(0,2*pi, sizeY), ...
                                       linspace(-pi/2, pi/2, sizeZ));
        
        scaleFactor = exp(-2*pi*sigma*(gridR./max(max(max(gridR)))));
        yLogPol = scaleFactor .* ifftnshift(X, 1:3);
        
        Y = convToCart(yLogPol, gridR, gridPhi, gridTheta, gridCartX, gridCartY, gridCartZ, interpMethod);
        
    case 'iD-AFMT'    
        
    otherwise 
    error('Unknown type of algorithm!')
    
end


    %% helper functions:
    function [xLogPol, gridR, gridPhi, gridTheta] = convToLogPol( xCart, sizeR, sizePhi, sizeTheta, interpMethod )
        [m, n, p]  = size(xCart);
        
        %Compute center positions oM and oN:
        oM = fix(m/2)+1;
        oN = fix(n/2)+1;
        oP = fix(p/2)+1;
        
        %Cartesian sample points:
        [gridX, gridY, gridZ] = meshgrid( -oN+1:n-oN, fliplr(-oM+1:m-oM), -oP+1:p-oP );
        
        %Log-polar plane sample points:
        maxR = sqrt(max(abs(oN-1), abs(n-oN))^2 + max(abs(oM-1), abs(m-oM))^2 + max(abs(oP-1), abs(p-oP))^2);
        [gridR, gridPhi, gridTheta]   = meshgrid( log(maxR)/sizeR:log(maxR)/sizeR:log(maxR), ...
                                      linspace(0,2*pi,sizePhi), ...
                                      linspace(-pi/2, pi/2, sizeTheta));
        
        %Perform conversion to log-polar coordinates, using linear interpolation:
        convGridX   = exp(gridR) .* cos(gridPhi) .* cos(gridTheta);
        convGridY   = exp(gridR) .* sin(gridPhi) .* cos(gridTheta);
        convGridZ   = exp(gridR) .* sin(gridTheta);
        xLogPol     = interp3(gridX, gridY, gridZ, xCart, ...
                        convGridX, convGridY, convGridZ, interpMethod, extrapVal);
    end

    function xCart = convToCart( xLogPol, gridR, gridPhi, gridTheta, gridX, gridY, gridZ, interpMethod )        
        %Perform conversion to cartesian coordinates:
%         convGridX   = reshape(exp(gridR) .* cos(gridPhi) .* cos(gridTheta),[], 1);
%         convGridY   = reshape(exp(gridR) .* sin(gridPhi) .* cos(gridTheta),[], 1);
%         convGridZ   = reshape(exp(gridR) .* sin(gridTheta),[], 1);
%         F = TriScatteredInterp(convGridX, convGridY, convGridZ, reshape(xLogPol, [], 1), interpMethod);
%         xCart       = F(gridX, gridY, gridZ);
        
        %interpolation in r-phi plane (worse but fast):
        convGridR       = sqrt(abs(gridX).^2 + abs(gridY).^2 + abs(gridZ).^2);
        convGridPhi     = atan(gridY./gridX);
        convGridTheta   = asin(gridZ./convGridR);
        %unwrap:
        qA = logical(logical((gridX >= 0)) .* logical((gridY<=0)));
        qB = logical(logical((gridX < 0)) .* logical((gridY<=0)));
        qC = logical(logical((gridX < 0)) .* logical((gridY>0)));
        qD = logical(logical((gridX >= 0)) .* logical((gridY>0)));
        convGridPhi(qA) = abs(convGridPhi(qA));
        convGridPhi(qB) = pi-convGridPhi(qB);
        convGridPhi(qC) = pi+abs(convGridPhi(qC));
        convGridPhi(qD) = 2*pi-convGridPhi(qD);
        convGridPhi(isnan(convGridPhi)) = 0;
        convGridTheta(isnan(convGridTheta)) = 0;
        
        xCart       = interp3(exp(gridR), gridPhi, gridTheta, xLogPol, ...
                        convGridR, convGridPhi, convGridTheta, interpMethod, extrapVal);
        
        
    end

end

