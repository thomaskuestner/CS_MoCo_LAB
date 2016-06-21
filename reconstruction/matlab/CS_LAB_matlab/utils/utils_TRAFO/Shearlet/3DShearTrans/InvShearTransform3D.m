function xRec = InvShearTransform3D(shCoeff)
level=size(shCoeff.D,2);


%Do the Band Pass of noisy Data
partialBP=cell(1,level);
BP=cell(1,level+1);

  %Assuming different band has same size of coeff
%     BP{l}=zeros(size(shCoeff.D{1,1}{1,1}));
  for pyrCone=1:3
    partialBP{pyrCone}=ShRec(shCoeff.D(pyrCone,:));
  end
  
 for l=1:level      
  %% Assuming different pyramidal zone have same shCoeff size at different 
  %%level
  BP{l}=zeros(size(partialBP{1}{l}));
  for pyrCone =1:3
   BP{l}=BP{l}+ partialBP{pyrCone}{l};
  end
end 
  
BP{level+1}=shCoeff.A;
%Reconstruct
xRec=DoPyrRec(BP);