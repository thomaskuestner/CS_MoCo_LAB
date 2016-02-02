function num_cores( desiredCores,flag_parpool )
% for implicit parallel processing
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

availableCores = maxNumCompThreads;
if availableCores <= desiredCores
    maxNumCompThreads(availableCores);
else
    maxNumCompThreads(desiredCores);
end;

%for explicit parallel processing:
if flag_parpool
    c = parcluster();
    p = gcp('nocreate');
    if isempty(p)
        p = parpool(c,desiredCores-1); % -1 for work on local pc
    end;
end;

end

