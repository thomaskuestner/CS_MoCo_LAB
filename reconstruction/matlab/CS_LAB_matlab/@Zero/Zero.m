classdef Zero < CSMaster
    % wrapper class for zero-padded reconstruction
    %
    % (c) Thomas Kuestner 
    % ---------------------------------------------------------------------
    
    properties      
        
    end
        
    methods
        function obj = Zero(measPara,espresso)
            obj = obj@CSMaster(1, 1, [], [], measPara); 
            obj.espresso = espresso;
        end
    end
    
    
end

