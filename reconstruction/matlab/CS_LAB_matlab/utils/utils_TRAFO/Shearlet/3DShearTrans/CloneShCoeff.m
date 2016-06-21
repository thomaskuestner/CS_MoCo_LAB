function shCoeffNew = CloneShCoeff(shCoeff,value)
shCoeffNew.D=cell(size(shCoeff.D));
shCoeffNew.A=shCoeff.A *value;

for l=1:size(shCoeff.D,2)
 for c=1:size(shCoeff.D,1)
  [L2 L1]=size(shCoeff.D{c,l});
  shCoeffNew.D{c,l}=cell(L2,L1);  
  for l2 =1:L2
    for l1=1:L1
      shCoeffNew.D{c,l}{l2,l1}=ones(size(shCoeff.D{c,l}{l2,l1}))*value;
    end
  end
 end
end

