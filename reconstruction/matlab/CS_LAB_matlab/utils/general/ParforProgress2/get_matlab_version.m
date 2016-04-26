function matlab_version = get_matlab_version()

    % old, slow version.
    %     v = ver('Matlab');
    %     matlab_version = str2double(v.Version);
     
    v = sscanf(version, '%d.%d.%d');
    matlab_version = [1 0.01 0.001] * v;
    
    % R2008a -> 7.06
    % R2008b -> 7.07
    % R2009a -> 7.08
    % R2009b -> 7.09
    % R2010a -> 7.10
    % R2010b -> 7.11
    % R2011a -> 7.12
    % R2011b -> 7.13 ?
    
end
%% EOF
