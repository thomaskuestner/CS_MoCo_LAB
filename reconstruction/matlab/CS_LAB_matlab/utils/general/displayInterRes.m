function displayInterRes( imgIn )
%DISPLAYINTERRES display intermediate results
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

% set all parameters here inside

persistent paraIR;

if(ischar(imgIn))   
    clear paraIR;
    return;
end

if(~isfield(paraIR,'dim'))
    paraIR.dim = ndims(imgIn);
    if(paraIR.dim == 3)
        % y-x-cha
        paraIR.subs = {':',':',':'};
        paraIR.permRule = [];
        
    elseif(paraIR.dim == 4)
        % 3D: y-z-x-cha
        
        paraIR.subs = {':', ':', ceil(size(imgIn,2)/2), ':'};
        paraIR.permRule = [1 3 2 4]; % => y-x-z-cha
        
%         % OR
%         % 2Dt: t-y-x-cha
%         paraIR.subs = {':', ':', ceil(size(imgIn,1)/2), ':'};
%         paraIR.permRule = [2 3 1 4]; % => y-x-t-cha
        
    elseif(paraIR.dim == 5)
        % t-y-z-x-cha
        paraIR.subs = {':', ':', ceil(size(imgIn,3)/2), ceil(size(imgIn,1)/2), ':'};
        paraIR.permRule = [2 4 3 1 5]; % => y-x-z-t-cha
    end
    
    paraIR.range = [0 1];
    paraIR.iter = 1;
    paraIR.normF = [];
    paraIR.screensize = get(0,'ScreenSize');
    
    paraIR.plotNorm = false;
end

if(~isempty(paraIR.permRule))
    imgIn = permute(imgIn,paraIR.permRule);
end
imgIn = imgIn(paraIR.subs{:});
imgIn = ((imgIn - min(imgIn(:))) * (paraIR.range(2)-paraIR.range(1)))./(max(imgIn(:)) - min(imgIn(:)));

% show channel combined abs image
imgIn = sqrt(sum(abs(imgIn).^2, ndims(imgIn)));
paraIR.normF(end+1) = norm(imgIn(:));

paraIR.hfig = figure(999);
set(paraIR.hfig,'position',[paraIR.screensize(1)+5 paraIR.screensize(2)+50 paraIR.screensize(3)-10 paraIR.screensize(4)-130]);
if(paraIR.plotNorm)
    subplot(2,1,1);
end
imagesc( imgIn, paraIR.range );
title(sprintf('\\bfiteration %g \nFrobenius norm: %g', paraIR.iter, norm(imgIn(:))))
axis image
colormap(gray(256))
% drawnow
if(paraIR.plotNorm)
    subplot(2,1,2);
    plot(paraIR.normF);
    xlabel('Iteration');
    ylabel('Frobenius norm');
    title('Frobenius norm behaviour');
end
drawnow
pause(0.5);

paraIR.iter = paraIR.iter + 1;

end

