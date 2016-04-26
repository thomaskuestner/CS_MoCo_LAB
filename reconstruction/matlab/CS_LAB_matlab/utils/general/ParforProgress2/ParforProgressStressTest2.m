function ParforProgressStressTest2(N)
% Stress test for 'ParforProgressStarter2'. In case of timeouts, you can
% use this function to determine how many simultaneous connections your
% computer can handle, and adjust the ppm.increment() call accordingly.
%
% Copyright (c) 2010-2012, Andreas Kotowicz
%
%%
    if nargin < 1
        N = 10000;
    end
    
    do_debug = 1;
    
    %% initialize ParforProgress monitor
    
    try % Initialization
        ppm = ParforProgressStarter2('test task - ParforProgressStarter2', N, 0.1, do_debug);
    catch me % make sure "ParforProgressStarter2" didn't get moved to a different directory
        if strcmp(me.message, 'Undefined function or method ''ParforProgressStarter2'' for input arguments of type ''char''.')
            error('ParforProgressStarter2 not in path.');
        else
            % this should NEVER EVER happen.
            msg{1} = 'Unknown error while initializing "ParforProgressStarter2":';
            msg{2} = me.message;
            print_error_red(msg);
            % backup solution so that we can still continue.
            ppm.increment = nan(1, N);
        end
    end

    %% execute dummy loop - replace 'rand(1, 10000)' with your function.
    
    t0 = tic();
    
    parfor i = 1 : N
        rand(1, 10000);
        ppm.increment(i); %#ok<PFBNS>
    end
    
    total_time = toc(t0);    
    
    %% clean up
    try % use try / catch here, since delete(struct) will raise an error.
        delete(ppm);
    catch me %#ok<NASGU>
    end
    
    %% show runtime

    disp(['ParforProgressStressTest2 running time: ' num2str(total_time) 's.']);
    
end
%% EOF
