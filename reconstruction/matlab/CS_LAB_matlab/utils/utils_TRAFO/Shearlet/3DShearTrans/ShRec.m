%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L= ShRec(shCoeff)
% Generates bandpass data as per the passed shearlet coefficient
%Input:   
%        shCoeff    : Shearlet Coefficient as obtained in Rec phase
%Output: 
%        L          : Reconstructed BandPass Data . Reconstuction
%                   : Can be Partial as per passed shearlet coeifficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L= ShRec(shCoeff)
level = length(shCoeff);
L=cell(1,level);

for l= 1:level
   [sizel2,sizel1] =size(shCoeff{l});
   D=zeros(size(shCoeff{l}{1,1}));
   for  l1= 1:sizel1
     for l2 =1:sizel2
      D=D+shCoeff{l}{l2,l1};
     end
   end 
   L{l}=D;
end