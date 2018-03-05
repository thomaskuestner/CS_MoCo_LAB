function num_cores_del( flag_parpool )
% deletes pool, if one was created
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

if flag_parpool
    delete(p);
else
    p = gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end;
end;

end

