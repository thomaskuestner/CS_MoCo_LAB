function [winout] = windowND(fhandle,n,winopt,para)
%WINDOWND create a N-dimensional window
% input:
% fhandle       window type
% n             1x1, 1x2 or 1x3 window width
% winopt        'symmetric' or 'periodic'
% para          additional window parameter

if(nargin < 3)
    winopt = 'symmetric';
end

if(nargin < 4)
    para = [];
end

if(size(n,2) > 3)
    error('windowND(): just 1, 2 or 3 dimensional window(s) allowed');
end

% if(ischar(fhandle))
% %     strcmp(fhandle,'@lanczos') || strcmp(fhandle,'@welch') || strcmp(fhandle,'@cosine') || strcmp(fhandle,'@planckTaper'))
%     w = windowOwn(fhandle,n,winopt,para);
% else

fname = functions(fhandle);
fname = fname.function;

w = cell(1,2);
for i=1:length(n)
    
    switch fname
        case {'rectwin','barthannwin','bartlett','bohmanwin','parzenwin','triang'}
            w{i} = window(fhandle,n(i));
        case {'blackman','blackmanharris','hamming','hann','flattopwin','nuttallwin'}
            w{i} = window(fhandle,n(i),winopt);
        case {'chebwin','gausswin','kaiser','tukeywin'}
            w{i} = window(fhandle,n(i),para);
        case {'lanczos','welch','cosine','planckTaper','kaiserbessel'}
            w{i} = windowOwn(fname,n(i),winopt,para);
    end
end

%     if(nargin < 3)
%         w = window(fhandle,n);
%     else
%         w = window(fhandle,n,winopt);
%     end
% end

if(size(n,2) == 3)
    w2D = w{1}(:) * w{3}(:).';
    winout = zeros(n(1), n(2), n(3));
    for i=1:n(2)
       winout(:, i, :) = w2D(:,:) * w{2}(i); 
    end     
elseif(size(n,2) == 2)
    winout = w{1}(:) * w{2}(:).';
else
    winout = w{1}(:); % * w{1}(:).';
end

end

function w = windowOwn(fhandle,n,symm,para)

% if(nargin < 3)
%     symm = 'symmetric';
% else
%     symm = winopt;
% end

% if(nargin < 4)
%     para = [];
% end

if(strcmp(symm,'symmetric'))
    L = n-1;
elseif(strcmp(symm,'periodic'))
    L = n;
else
    error('windownOwn(): Unknown symmetry type');
end


w = (0:n-1); % - L/2;
if(strcmp(fhandle,'lanczos'))
    w = sinc(2*w(:)./(n-1) - 1);
    
elseif(strcmp(fhandle,'welch'))
    w = 1 - ((w(:) - (n-1)/2)/((n+1)/2)).^2;

elseif(strcmp(fhandle,'cosine'))
    w = sin((pi*w(:))./(n-1));
    
elseif(strcmp(fhandle,'planckTaper'))
    if(isempty(para))
        error('windowOwn(): Please specifiy epsilon parameter for Planck-Taper window');
    end
    tl = max([floor(para*(n-1)),2]); %+1 ;
    th = min([(n-1)-tl+2,length(w)-1]); % + 1;
        
    Zp = (n-1)*para*(1./w(2:tl) + 1./(w(2:tl)-(n-1)*para));
    Zm = -(n-1)*para*(1./(w(th:end-1)-(n-1)+(n-1)*para) + 1./(w(th:end-1)-(n-1)));
    
%     Zp = 2*para*(1./(1+(2*t(1:tl)./(n-1))) + 1./(1-2*para+(2*t(1:tl)./(n-1))));
%     Zm = 2*para*(1./(1-(2*t(th:end)./(n-1))) + 1./(1-2*para-(2*t(th:end)./(n-1))));
    
    w = [1./(exp(Zp)+1), ones(1,th-tl+1), 1./(exp(Zm)+1)]';
    
    
%     bounds = [0, round(para*(n-1)), round((1-para)*(n-1)), n-1]+1;
%     Zp = 2*para*(1./(1+2*w(bounds(1):bounds(2))./(n-1)) + 1./(1-2*para+2*w(bounds(1):bounds(2))./(n-1)));
%     Zm = 2*para*(1./(1-2*w(bounds(3)+1:bounds(4))./(n-1)) + 1./(1-2*para-2*w(bounds(3)+1:bounds(4))./(n-1)));
%     
%     w = zeros(1,n);
%     w(bounds(1):bounds(2)) = 1./(exp(Zp)+1);
%     w(bounds(2)+1:bounds(3)) = 1;
%     w(bounds(3)+1:bounds(4)) = 1./(exp(Zm)+1);
    
elseif(strcmp(fhandle,'kaiserbessel'))
    if(isempty(para))
        error('windowOwn(): Please specifiy beta parameter for Kaiser-Bessel window');
    end
    
    tl = floor(n/2);
    
    wK = window(@kaiser,tl,para);
    
    sumwK = sum(wK(:));
    cumsumwK = cumsum(wK(:));
    
    w = [sqrt(cumsumwK(1:tl)./sumwK); sqrt(cumsumwK(tl:-1:1)./sumwK)];
    
else
    error('windowOwn(): Undetermined window type');
end

end