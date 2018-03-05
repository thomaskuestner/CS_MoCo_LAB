function dispProgress(name, value, maxvalue)
%DISPPROGRESS display calculation progress
% initialize and update progress bar
%
% (c) Thomas Kuestner
% ---------------------------------------------------------------------

persistent prop

if(nargin < 3)
    maxvalue = 0;
end

if(~exist('prop','var') || isempty(prop) || strcmp(name,'InitDispProgress'))
    if(isstruct(value))
        prop = value;
        return;
    else
%         error('dispProgress(): Unknown data structure for property'); 
        return; % shortcut for GA simulation
    end
end

if(prop.flagParallel)
    if(verLessThan('matlab','8.4'))
        if(matlabpool('size') > 0)
            return; % shortcut for parallel computing
        end
    else
        if(~isempty(gcp('nocreate')))
            return % shortcut for parallel computing
        end
    end
end

% if(matlabpool('size') > 0)
%     openPool = true;
%     currpath = fileparts(fileparts(mfilename('fullpath')));
%     a = load([currpath, filesep, 'prop.mat']);
% else
%     openPool = false;
% end

if(prop.flagDisp)
    % GUI

    if(~prop.openPool)
        if(prop.disp.openBar(strcmp(prop.disp.names,name)))
            % progress bar already exists
            if(ischar(value) && strcmp(value,'Close'))
                multiWaitbar(name,value);
                prop.disp.openBar(strcmp(prop.disp.names,name)) = false;
            else
                cval = round(value*100);
                if(cval == 0), cval = max(ceil(value*100),1); end;
                if(cval > 100), cval = 100; end;
                abort = multiWaitbar(name, value, 'Color', prop.disp.colors(cval,:));
                if(abort)
                    multiWaitbar('CLOSEALL');
                    error('User Termination of CS reconstruction!');
                end     
            end
        else
            % progress bar needs initialization
            if(maxvalue > 1)
                multiWaitbar(name,value,'CancelFcn', @(a,b) disp( ['Cancel ',a] ));
                prop.disp.openBar(strcmp(prop.disp.names,name)) = true;
            end
        end

    else
        % parallel computation
        if(prop.disp.openBar(strcmp(prop.disp.names,name)))
            % progress bar already exists
            if(ischar(value) && strcmp(value,'Close'))
                try
                    delete(prop.disp.ppm{strcmp(prop.disp.names,name)});
                catch err
                end
                prop.disp.openBar(strcmp(prop.disp.names,name)) = false;
            else
                prop.disp.ppm{strcmp(prop.disp.names,name)}.increment(value);
            end
        else
            % progress bar needs initialization
            if(maxvalue > 1)
                helper = ParforProgressStarter2(name, maxvalue, 1);
                prop.disp.ppm{strcmp(prop.disp.names,name)} = helper;
                prop.disp.openBar(strcmp(prop.disp.names,name)) = true;
            end
        end
%         save([currpath, filesep, 'prop.mat'], '-struct', 'a');
    end
    
else
    % console
    
    if(ischar(value) && strcmp(value,'Close'))
        if(any(strcmp(name,{'CG','SENSE - Time','SENSE - Kernel','SENSE - Concatenate','Extracting K-Space','GRAPPA calibration','Partial Fourier','ESPReSSo','Trafo','Log barrier','Newton','Proximal Average', 'ADMM'})))
            if(prop.openPool)
                fprintf([repmat('\b',1,length(name)+1)]);
            else
                fprintf([repmat('\b',1,length(prop.lastOut))]);
            end
        end
    
    elseif(value == 0)
        % init step
        prop.maxvals(strcmp(prop.disp.names,name)) = maxvalue;
        if(strcmp(name,'FOCUSS'))
%             prop.lastOut = sprintf('%s %3d%%',name,value);
            prop.lastVal = 0;
            fprintf('%s [%s]\n%s [', name, repmat('- ',1,25),repmat(' ',1,length(name)));
        elseif(strcmp(name,'CG'))
%             prop.lastOut = sprintf('%s\t\t%s %3d%%\n',prop.lastOut,name,value);
%             if(~prop.openPool)
%                 fprintf('%s',prop.lastOut);
%             else
%                 prop.lastOut = '';
%             end
            prop.lastVal = 0;
            fprintf('%s [%s]\n%s [', name, repmat('- ',1,25),repmat(' ',1,length(name)));
        elseif(any(strcmp(name,{'Extracting K-Space','GRAPPA calibration','Partial Fourier','ESPReSSo','Trafo','Log barrier','Newton','Proximal Average', 'ADMM'})))
%             if(prop.openPool)
%                 prop.lastOut = sprintf('%s\n',name);
%             else
%                 prop.lastOut = sprintf('%s %3d%%\n',name,value);
%             end
%             fprintf('%s',prop.lastOut);
            prop.lastVal = 0;
            fprintf('%s [%s]\n%s [', name, repmat('- ',1,25),repmat(' ',1,length(name)));
        end
        
    else
        if(~prop.openPool)
            switch name
                case {'Slices','Channels','Repetitions','POCS','Time','Frequency','Averages'}
                    outstring = sprintf('%s %d/%d\n', name, value*prop.maxvals(strcmp(prop.disp.names,name)), prop.maxvals(strcmp(prop.disp.names,name)));
                    delString = 0;

                case {'Extracting K-Space','GRAPPA calibration','Partial Fourier','ESPReSSo','Trafo','Log barrier','Newton','Line Search','Proximal Average', 'ADMM'}
%                     outstring = sprintf('%s %3d%%\n', name, value*100); 
%                     delString = length(prop.lastOut);

                    value = round(value*100);
                    if(mod(value,2) ~= 0)
                        value = value-1;
                    end
                    value = value/2;
                    
                    if(value >= 50)
                        outstring = sprintf('*]\n');
                    else
                        if(value == prop.lastVal)
                            outstring = '';
                        else
                            outstring = '*'; 
                        end
                    end
                    prop.lastVal = value;
                    delString = 0;          
                case {'FOCUSS', 'CG'}

%                     if(strcmp(name,'FOCUSS'))
%                         outstring = sprintf('%s %3.0f%%%s',name,value*100,prop.lastOut(12:end));
%                     else
%                         outstring = sprintf('%s\t\t%s %3.0f%%\n', prop.lastOut(1:11),name,value*100);
%                     end
% 
%                     if(~isempty(prop.lastOut))
%                         delString = length(prop.lastOut);
%                     else
%                         delString = 0;
%                     end
                    value = round(value*100);
                    if(mod(value,2) ~= 0)
                        value = value-1;
                    end
                    value = value/2;
                    
                    if(value >= 50)
                        outstring = sprintf('*]\n');
                    else
                        if(value == prop.lastVal)
                            outstring = '';
                        else
                            outstring = '*'; 
                        end
                    end
                    prop.lastVal = value;
                    delString = 0;                    

                case 'SENSE - Slice'
                    outstring = sprintf('%s %d/%d\n', name, value*prop.maxvals(strcmp(prop.disp.names,name)), prop.maxvals(strcmp(prop.disp.names,name)));
                    delString = 0;

                case 'SENSE - Time'
                    outstring = sprintf('%s %d/%d\n', name, value*prop.maxvals(strcmp(prop.disp.names,name)), prop.maxvals(strcmp(prop.disp.names,name)));

                    if(~isempty(prop.lastOut))
                        delString = length(prop.lastOut);
                    else
                        delString = 0;
                    end

                case {'SENSE - Kernel', 'SENSE - Concatenate'}
                    outstring = sprintf('%s %3.0f%%\n',name,value);
                    if(~isempty(prop.lastOut))
                        delString = length(prop.lastOut);
                    else
                        delString = 0;
                    end

                otherwise
                    error('dispProgress(): Undetermined display name');              
            end

            prop.lastOut = outstring;

            fprintf([repmat('\b',1,delString),'%s'],outstring);
        end
    end
end

end

