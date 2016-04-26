classdef ParforProgressConsole2 < handle
% do NOT use this object by itself!
% use ParforProgressStarter() instead.    
% not using saveobj and loadobj, because old matlab versions have problems
% with it.
% Copyright (c) 2010-2012, Andreas Kotowicz

    properties (GetAccess = private, SetAccess = private)
        message
        nbr_files
        fraction
        start_time
    end

    methods
        
        function o = ParforProgressConsole2(s, n, percentage)
            
            if nargin < 3
                percentage = 0.05;
            end
            
            if nargin < 2
                n = 5;
            end
            
            if nargin < 1
                s = 'test';
            end
            
            o.fraction = round(n * percentage);
            o.nbr_files = num2str(n);
            o.message = s;
            o.start_time = tic;

        end

        function increment(o, i)
            if mod(i, o.fraction) == 0
                disp([num2str(i) ' / ' o.nbr_files ' ' o.message '.']);
            end
        end
        
        function delete(o)
            % if you have muliple workers, then this message will show up
            % for each object.
            if ~isempty(o.start_time)
                blub = toc(o.start_time);
                disp(' ');
                disp(['  >> execution time was ' num2str(blub) 's.']);
                o.start_time = [];
            end
        end

    end
    
end
%% EOF
