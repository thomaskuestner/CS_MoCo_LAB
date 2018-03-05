function F= GetCharacteristicBasedFilter(level,dBand,dataClass)

F=cell(3,level);  
for l=1:level  
  cubeSize=8;
  [P ,PF]=GeneratePyramidSection(cubeSize);
  
  numDir=dBand{level}{l}(1,1);    
  shift=cubeSize/numDir;
  F{1,l}=cell(numDir,numDir);
  F{2,l}=cell(numDir,numDir);
  F{3,l}=cell(numDir,numDir);
  SF=PF{1}+PF{2}+PF{3};
  mRadial =ones(shift,cubeSize,shift);
  
  A=zeros(cubeSize,cubeSize,cubeSize);
  
   for l2=1:numDir
    for l1=1:numDir
      for c=1:3
      mRadialIdx=[(l1-1)*shift+1  l1*shift ;1 cubeSize ; (l2-1)*shift+1  l2*shift  ];
      
      F{c,l}{l2,l1}=PolarToRec(mRadialIdx,mRadial,cubeSize,P{c},SF);
       A=A+F{c,l}{l2,l1};
%       F{c,l}{l2,l1}=squeeze(real(fftshift(ifftn(fftshift(F{c,l}{l2,l1})))));

%       F{2,l}{l2,l1}=PolarToRec(mRadialIdx,mRadial,cubeSize,P{2},SF);
%       A=A+F{2,l}{l2,l1};
%       F{2,l}{l2,l1}=real(fftshift(ifftn(fftshift(F{2,l}{l2,l1}))));     
%       F{3,l}{l2,l1}=PolarToRec(mRadialIdx,mRadial,cubeSize,P{3},SF);
%       A=A+F{3,l}{l2,l1};
%       F{3,l}{l2,l1}=real(fftshift(ifftn(fftshift(F{3,l}{l2,l1}))));
      end
    end
   end
for c=1:3
  for l2=1:numDir
    for l1=1:numDir
       F{c,l}{l2,l1}=F{c,l}{l2,l1}./A;
        F{c,l}{l2,l1}=squeeze(real(fftshift(ifftn(fftshift(F{c,l}{l2,l1})))));   
     end
   end
 end
end
                                           

Temp=F(1,:);
F(1,:)=F(2,:);
F(2,:)=Temp;
