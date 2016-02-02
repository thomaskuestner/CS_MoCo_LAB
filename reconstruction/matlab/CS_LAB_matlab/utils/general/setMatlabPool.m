function setMatlabPool(val, prop)
% change status of matlab pool
% val: open|close
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

if(nargin == 0)
    val = 'close';
end

if(prop.flagParallel)
    % parallel computation toolbox exists
    if(strcmp(val,'open'))
        if(~prop.openPool)
            if(matlabpool('size') > 0)
                warning('Closing an unopened MatLab pool');
                matlabpool close
            end
            matlabpool open
            prop.openPool = true;
        else
            warning('A currently used MatLab pool is already open');
        end


    elseif(strcmp(val,'close'))
        if(prop.openPool)
            matlabpool close
        else
            if(matlabpool('size') > 0)
                warning('Closing an unopened MatLab pool');
                matlabpool close
            end
        end
        prop.openPool = false;
    else
        error('setMatlabPool(): Undetermined value');
    end
end

end
