function [x, active] = pool_size()
%POOL_SIZE - temporary hack to return size of current MATLABPOOL
%
% see:
% http://www.mathworks.com/support/solutions/en/data/1-5UDHQP/index.html?product=DM&solution=1-5UDHQP

    %% get matlabpool size - this should work for ALL matlab versions
    try
        session = com.mathworks.toolbox.distcomp.pmode.SessionFactory.getCurrentSession;

        if ~isempty( session ) && session.isSessionRunning() && session.isPoolManagerSession()
            client = distcomp.getInteractiveObject();
            if strcmp( client.CurrentInteractiveType, 'matlabpool' )
                x = session.getLabs().getNumLabs();
            else
                x = 0;
            end
        else
            x = 0;
        end

    catch me %#ok<NASGU>
        % standalone matlab installations might not have the appropriate
        % com module installed.
        x = 0;
    end
    

    %% matlabpool active?
    active = 0;
    if x > 0
        active = 1;
    end

end
%% EOF
