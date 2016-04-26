% This  function adds to 3D shearlet    coefficients as per weight Assuming
% same number of level and directional band, coarase level coefficient size
%
%TODO Add checks for size and direction and level mismatch
function shCoeff1 = AddShcoeff(shCoeff1,shCoeff2,weight1,weight2)

shCoeff1.A=  weight1*shCoeff1.A+ weight2*shCoeff2.A;
for l=1:size(shCoeff1.D,2)
 for c=1:size(shCoeff1.D,1)
  [L2 L1]=size(shCoeff1.D{c,l});
%   shCoeffNew.D{c,l}=cell(L2,L1);
  for l2 =1:L2
    for l1=1:L1
      shCoeff1.D{c,l}{l2,l1}= weight1*shCoeff1.D{c,l}{l2,l1}+ ... 
         weight2*shCoeff2.D{c,l}{l2,l1};
    end
  end
 end
end
