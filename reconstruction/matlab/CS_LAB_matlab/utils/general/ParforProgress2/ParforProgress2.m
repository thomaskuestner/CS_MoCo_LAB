%
% do NOT use this object by itself!
% use ParforProgressStarter2() instead.
% 
% most of the code here is from:
% http://www.mathworks.com/matlabcentral/fileexchange/24594-parfor-progress-monitor
%
% ParforProgress2 - M object to make 'ParforProgressClient2' and 
% 'ParforProgressServer2' objects easier to use. Create one of these on the 
% client outside your PARFOR loop with a name for the window. Pass the 
% object to the PARFOR loop, and have the workers call "increment" at the 
% end of each iteration. This sends a notification back to the server which 
% then updates the GUI.
%
% Example:
%
% N = 100;
% ppm = ParforProgress2('my task', N);
% parfor i = 1 : N
%     rand(1);
%     ppm.increment(i);
% end
% delete(ppm);
%
%
% Copyright (c) 2010-2012, Andreas Kotowicz

classdef ParforProgress2 < handle

    properties (GetAccess = private, SetAccess = private)
        Port
        HostName
        OldVersion = 0;
        DEBUG = 0;
    end
    
    properties (Transient, GetAccess = private, SetAccess = private)
        JavaBit
    end
    
    methods (Static)
        function o = loadobj(X)
        % Once we've been loaded, we need to reconstruct ourselves
        % correctly as a client
            o = ParforProgress2( {X.HostName, X.Port, X.OldVersion, X.DEBUG} );
            % disp('loading');
        end
    end
    
    methods
        function o = ParforProgress2(s, n, percentage, do_debug, use_gui)
        % ParforProgress Build a Parfor Progress Monitor
        % Use the syntax: ParforProgress( 'Window Title', N, percentage, do_debug, use_gui )
        % where N is the number of iterations in the PARFOR loop
        
            % initalize client
            if nargin == 1 && iscell(s)
                % "Private" constructor used for the clients
                o.JavaBit   = ParforProgressClient2.createClient(s{1}, s{2}, s{4});
                o.Port      = [];
             
            % initialize server
            elseif (nargin == 5 || nargin == 4 || nargin == 3 || nargin == 2)
                
                if nargin < 5
                    use_gui = 1;
                end
                
                if nargin < 4
                    do_debug = 0;
                end
                
                if nargin < 3
                    percentage = 0.1;
                end
                
                % check for old matlab versions.
                if get_matlab_version < 7.07
                    o.OldVersion = 1;
                end

                o.JavaBit   = ParforProgressServer2.createServer(s, n, percentage, use_gui);
                o.Port      = double(o.JavaBit.getPort());
                
                % Get the client host name from pctconfig - needs
                % distcomp toolbox.
                % cfg         = pctconfig;
                % o.HostName  = cfg.hostname;
                
                % gethostname() is also problematic, because laptop might
                % not be connect to network (will have local IP only).
                %o.HostName = gethostname();
                
                address = java.net.InetAddress.getLocalHost;
                o.HostName = char(address.getHostAddress);

                o.DEBUG = do_debug;
                
            else
                error( 'Public constructor is: ParforProgress2(''Text'', N, percentage, do_debug, use_gui)' );
            end
        end
        
        function X = saveobj(o)
        % keep port, hostname, matlab version and debug flag
            X.Port     = o.Port;
            X.HostName = o.HostName;
            X.OldVersion = o.OldVersion;
            X.DEBUG = o.DEBUG;
            % disp('saving');
        end
        
        function increment(o, i) %#ok<INUSD>
        % i is a fake input so we stay compatible with
        % "ParforProgressConsole2.m"
        
            % Something is wrong with matlab 2008a and the saveobj /
            % loadobj methods. There is no 'JavaBit' in the 'o' if you
            % have matlabpool enabled.
            if o.OldVersion == 0
                o.JavaBit.increment();
            end
            
        end
        
        function delete(o)
        % Close the UI
            % somethings wrong with matlab 2008a and the saveobj /
            % loadobj methods. there is no 'JavaBit' in the 'o' if you
            % have matlabpool enabled.
            if o.OldVersion == 0
                o.JavaBit.done();
            end
        end
        
    end
    
end
%% EOF
