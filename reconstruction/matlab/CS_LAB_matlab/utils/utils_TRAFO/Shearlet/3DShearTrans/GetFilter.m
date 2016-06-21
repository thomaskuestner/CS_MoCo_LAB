function F=GetFilter(filterName,level,dBand,filterSize,filterDilationType ,dataClass)

switch filterName    
    case 'meyer'           
          
        switch filterDilationType
              case '422'
                  F= GetMeyerBasedFilter(level,dBand, filterSize ,dataClass);
              otherwise
                  error('filter %s with dilation type %s not implemented',filterName,filterDilationType);
        end
        
        
    otherwise
         error('filter %s not implemented',filterName );
end
