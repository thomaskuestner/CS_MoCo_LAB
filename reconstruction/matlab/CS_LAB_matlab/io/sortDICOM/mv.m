function mv(arg1, arg2, options)
% mv(arg1, arg2, options)
% move arg1 to arg2
% if options == 'try', move can fail without throwing an error
if ~exist('options', 'var'); options = {}; end
if ~iscell(options); options = {options}; end
if (strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64'))
    arg1 = regexprep(arg1,'/','\');
    arg1 = regexprep(arg1, '\\{2,}', '\\');
    arg2 = regexprep(arg2,'/','\');
    arg2 = regexprep(arg2, '\\{2,}', '\\');
    %system(['move "' arg1 '" "' arg2 '"']);
    %arg1 = makeUNCpath(arg1);
    %arg2 = makeUNCpath(arg2);
    if ismember({'try'}, options)
        status = movefile(arg1, arg2);
    else
%         system(['move "' arg1 '" "' arg2 '" >nul']); % does e.g. not work with UNC pathnames
        movefile(arg1, arg2);
    end
else
    system(['mv ' arg1 ' ' arg2]);
end


function UNCpath = makeUNCpath(path)
% make UNC path to allow for long path names
if isempty(strfind(strrep(path, '/', '\'),'\'))
    UNCpath = ['\\?\' pwd '/' path];
elseif (length(path) > 5) && ~isequal(path(1:4), '\\?\')
    if isequal(path(2), ':')
        UNCpath = ['\\?\' path];
    else
        UNCpath = path;
        % TODO: compute missing part of drive letter / path if relative path
        % given
    end
else
    UNCpath = path;
    % TODO: compute missing part of drive letter / path if relative path
    % given
end
