function res = fftnshift(x,fftdim,scrambledim)
% FFT
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

    if(nargin == 1)
        % full N-D FFT
        fftdim = 1:ndims(x);
    end

%     if(exist(['utils',filesep,'general',filesep,'fftw_wisdom.mat'],'file') || exist([fileparts(mfilename('fullpath')),filesep,'fftw_wisdom.mat'],'file'))
%         load([fileparts(mfilename('fullpath')),filesep,'fftw_wisdom.mat']);
    if(any(fftdim > ndims(x)))
        fftdim = fftdim(fftdim <= ndims(x));
    end
    if(isempty(fftdim))
        res = x;
        return;
    end
    if(evalin('base','exist(''fftw_wisdom'',''var'')'))
        sizesCalculated = evalin('base', 'fftw_wisdom.sizesCalculated');
        wisdom_str = evalin('base', 'fftw_wisdom.wisdom_str');
        depth = sum(fftdim);
        if(length(fftdim) > 1), depth = depth + 1; end;
        if(depth <= size(sizesCalculated,3))
            sizesCalculated = sizesCalculated(:,fftdim,depth);
            xSize = size(x);
            line = ismember(sizesCalculated,xSize(fftdim),'rows');
            if(~isempty(line))
                fftw('wisdom',wisdom_str{line,depth});
            end
        end
    end    
    
    res = x;
    if(nargin < 3)
        % fftdim == scrambledim (scramble all directions)
        for i=1:length(fftdim)
%             res = 1/(size(res,fftdim(i))) * fftshift(fft(ifftshift(res,fftdim(i)),[],fftdim(i)),fftdim(i));
            res = fftshift(fft(ifftshift(res,fftdim(i)),[],fftdim(i)),fftdim(i));
        end
    else
        % partial FFT (just scramble specific directions)
        dim = unique([fftdim, scrambledim]);
        
        for i=1:length(dim)
            if(ismember(dim(i),fftdim) && ismember(dim(i),scrambledim))
%                 res = 1/(size(res,dim(i))) * fftshift(fft(ifftshift(res,dim(i)),[],dim(i)),dim(i));
                res = fftshift(fft(ifftshift(res,dim(i)),[],dim(i)),dim(i));
            elseif(ismember(dim(i),fftdim) && ~ismember(dim(i),scrambledim))
%                 res = 1/(size(res,dim(i))) * fft(res,[],dim(i));
                res = fft(res,[],dim(i));
            elseif(~ismember(dim(i),fftdim) && ismember(dim(i),scrambledim))
                % void fftshift(ifftshift())
            else
                error('ifftnshift(): Impossible logical expression');
            end
               
        end
    end
end


