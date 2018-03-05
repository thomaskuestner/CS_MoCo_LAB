%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function shCoeff= ShDec(pyrCone,F,L,dataClass)
% Computes shearlet coefficients
%Input:   
%        pyrCone    : pyramidal cone index (1 ,2 ,3)
%        F          : Windowing filter cell array
%        L          : BandPass data/Other preproccessed data
%        dataClass  : 'single' or 'double' 
%Output: 
%        shCoeff    : MultiScale Shearlt Coefficient
%                   : First level cell index gives scale starting from
%                   : from finest  level. Seond level cell index gives
%                   : coefficient as per index location in frequency domain
%                   : like shCoeff{1}{2,2} gives shearlet Coefficient
%                   : at finest level and band orientation L1=2,L2=2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shCoeff= ShDec(pyrCone,F,L,level,dataClass)
%level=length(L)-1;
shCoeff=cell(1,level);
for l=1:level  
 [sizel2,sizel1] =size(F{pyrCone,l});
    for l2=1:sizel2
      for l1=1:sizel1  
        if strcmp(dataClass,'single')
	       shCoeff{l}{l2,l1}=single(convnfft(L{l},F{pyrCone,l}{l2,l1},'same'));
        elseif strcmp(dataClass,'double')
          shCoeff{l}{l2,l1}=double(convnfft(L{l},F{pyrCone,l}{l2,l1},'same'));
        elseif strcmp(dataClass,'uint8')
          shCoeff{l}{l2,l1}=uint8(convnfft(L{l},F{pyrCone,l}{l2,l1},'same'));
        end
      end    
    end
end
