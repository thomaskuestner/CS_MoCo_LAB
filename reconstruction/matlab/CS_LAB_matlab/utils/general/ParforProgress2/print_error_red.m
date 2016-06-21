function print_error_red(msg, tab_before)
% show error in red

    if nargin < 1
        disp('usage: print_error_red(msg)');
        return;
    end
    
    if nargin < 2
        tab_before = 1;
    end

    if ~iscell(msg)
        disp('error message must be a cell array of strings.');
        return;
    end
    
    % make sure we have a column vector
    if size(msg, 1) > size(msg, 2)
        msg = msg';
    end

    %% create a fake error output
    if tab_before == 1
        for i = 1 : numel(msg)
            fprintf(2, '\t %s \n', msg{1, i});
        end
    else
        for i = 1 : numel(msg)
            fprintf(2, ' %s \n', msg{1, i});
        end        
    end
    
end
%% EOF
