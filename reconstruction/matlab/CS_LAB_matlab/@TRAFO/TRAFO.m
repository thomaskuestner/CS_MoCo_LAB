classdef TRAFO < handle
    % TRAFO Transformation class to overload mathematical operations
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------

    properties
        data;
        meta;
    end
    
    methods
        function obj = TRAFO(data,meta,paraTrafo,initialize)
            if(nargin < 3)
                obj.data = data;
                obj.meta = meta;
            else
                if(nargin < 4)
                    initialize = 'default';
                end
                [obj.data, obj.meta] = initTRAFO(data,paraTrafo,meta,initialize);
%                 if(isempty(obj.meta))
%                     obj.meta = meta;
%                 end
            end
        end
        
        function data = getData(obj)
            data = obj.data;
        end
        
        function meta = getMeta(obj)
            meta = obj.meta;
        end
           
    end
    
end

