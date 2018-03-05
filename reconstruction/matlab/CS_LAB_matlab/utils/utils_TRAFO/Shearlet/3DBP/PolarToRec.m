%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function mBlock= PolarToRec(mRadialSlice,cubeSide,cPyramid)
%% Generates the XYZ coordiante of Polar Grid in Cartesian Grid
%%Input:   mRadialSlice: 3-d Radial slcie matrix
%%        cubeSide: order of block matrix in 3-d to be assembled
%%        cPyramid: cell containing one of 3 pyramid information
%%Output: 
%%        mBlock: block matrix in 3-d to be assembled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mBlock= PolarToRec(mRadialIdx,mRadial,cubeSide,pyramid,F)
mBlock=zeros(cubeSide,cubeSide,cubeSide); 
for j=mRadialIdx(3,1):mRadialIdx(3,2)
  for i=mRadialIdx(1,1):mRadialIdx(1,2)
    for k=mRadialIdx(2,1):mRadialIdx(2,2)
      mBlock(pyramid.X(i,k,j),pyramid.Y(i,k,j),pyramid.Z(i,k,j))=...
       +mRadial(i-mRadialIdx(1,1)+1, k-mRadialIdx(2,1)+1, j-mRadialIdx(3,1)+1);
    end
  end
end    

% FToggle=F==0;
% F=F+FToggle;
%   mBlock=mBlock./F;



