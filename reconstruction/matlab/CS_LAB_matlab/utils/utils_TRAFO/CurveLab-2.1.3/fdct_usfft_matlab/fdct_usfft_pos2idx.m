function [a,b] = fdct_usfft_pos2idx(NX,NY,SX,SY,s,w,x,y)

%fdct_usfft_pos2idx.m - For a fixed scale and fixed direction, returns
%		the curvelet which is closest to a certain point on the image
%
% Inputs
%   NX,NY,SX,SY     Values returned by fdct_usfft_param
%   s               scale index
%   w               wedge (angular) index
%   x,y             position in image
%
% Outputs
%   a,b             Index of the curvelet at scale s and angle w which is nearest to (x,y)
%
  
  nx = NX{s}{w};  ny = NY{s}{w};
  bx = SX{s}{w}(1,1);  by = SY{s}{w}(1,1);
  sx = [SX{s}{w}(2,1)-bx, SY{s}{w}(2,1)-by];
  sy = [SX{s}{w}(1,2)-bx, SY{s}{w}(1,2)-by];
  tmp = 1 + [x-bx,y-by] / [sx; sy];
  a = round(tmp(1));  b = round(tmp(2));
  a = mod(a,nx);  b = mod(b,ny);
  if(a<1)    a = a+nx;  end
  if(a>nx)    a = a-nx;  end
  if(b<1)    b = b+ny;  end
  if(b>ny)    b = b-ny;  end
  
