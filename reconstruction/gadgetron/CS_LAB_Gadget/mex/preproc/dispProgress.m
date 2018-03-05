function dispProgress(name, value, ~)
%DISPPROGRESS display calculation progress
% initialize and update progress bar

persistent prop

if(nargin < 3)
    if(ischar(value) && strcmp(value,'Close'))
        close(prop.hw);
    else
        waitbar(value,prop.hw);
    end
else
    % init progress bar
    prop.hw = waitbar(value, name);
end

end

