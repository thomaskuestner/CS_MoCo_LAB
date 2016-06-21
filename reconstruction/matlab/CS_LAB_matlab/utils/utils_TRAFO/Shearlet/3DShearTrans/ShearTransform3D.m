function shCoeff = ShearTransform3D(X,shearingFilter)
level=size(shearingFilter,2);
shCoeff.D=cell(3,level);

%Do the Band Pass of noisy Data
BP=DoPyrDec(X,level); 
shCoeff.A=BP{1};
for pyrCone=1:3
  shCoeff.D(pyrCone,:) = ShDec(pyrCone,shearingFilter,BP,level,'single');
end
