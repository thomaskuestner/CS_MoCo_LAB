function [SX,SY,FX,FY,NX,NY] = fdct_usfft_param(C)

% fdct_usfft_param: gives the phase space location of each curvelet
%
% Input
%     C     Curvelet coefficients
%
% Output
%     SX    Cell array, SX{s}{w}(i,j), gives the row index of the center of the curvelet
%           index by (i,j)  in scale s and wedge w 
%     SY    Cell array, SY{s}{w}(i,j), gives the column index of the center of the curvelet
%           index by (i,j)  in scale s and wedge w 
%     FX    Cell array, FX{s}{w} gives the row index of the center frequency of the curvelet
%           in scale s and wedge w
%     FY    Cell array, FY{s}{w} gives the column index of the center frequency of the curvelet
%           in scale s and wedge w
%     NX    Cell array, NX{s}{w} gives the number of rows of the curvelet matrix
%           in scale s and wedge w
%     FY    Cell array, FY{s}{w} gives the number of columns of the curvelet matrix
%           in scale s and wedge w
  
  [m,n] = size(C{end}{1});
  nbs = length(C);
  L = 4;
  SX = cell(1,nbs);  SY = cell(1,nbs);
  FX = cell(1,nbs);  FY = cell(1,nbs);
  NX = cell(1,nbs);  NY = cell(1,nbs);
  
  %first
  SX{1} = cell(1,1);  SY{1} = cell(1,1);
  FX{1} = cell(1,1);  FY{1} = cell(1,1);
  NX{1} = cell(1,1);  NY{1} = cell(1,1);
  
  FX{1}{1} = 0;  FY{1}{1} = 0;
  [nx,ny] = size(C{1}{1});
  NX{1}{1} = nx;  NY{1}{1} = ny;
  sx = m/nx;  sy = n/ny;
  [TX,TY] = ndgrid(0:nx-1,0:ny-1);
  SX{1}{1} = 1 + TX*sx;  SY{1}{1} = 1 + TY*sy;
  
  ctr = m/(3*2^(nbs-2)); %since m==n here
  %middle ones
  for s=2:nbs-1
    nw = length(C{s});
    SX{s} = cell(1,nw);  SY{s} = cell(1,nw);
    FX{s} = cell(1,nw);  FY{s} = cell(1,nw);
    NX{s} = cell(1,nw);  NY{s} = cell(1,nw);
    cnt = 1;
    nwpq = nw/4;
    %quadrant 1
    for w=1:nwpq
      FX{s}{cnt} = -ctr;      FY{s}{cnt} = -ctr + (w-0.5)*2*ctr/nwpq;
      s2 = -FY{s}{cnt}/ctr;
      [nx,ny] = size(C{s}{cnt});
      NX{s}{cnt} = nx;      NY{s}{cnt} = ny;
      sx = [1/nx, 0];      sy = [-s2*1/ny, 1/ny];
      cx = nx/2+1;      cy = ny/2+1;
      [IX, IY] = ndgrid(1:nx, 1:ny);
      SX{s}{cnt} = 1 + m * mod(sx(1)*(IX-cx)+sy(1)*(IY-cy)+0.5,1);
      SY{s}{cnt} = 1 + n * mod(sx(2)*(IX-cx)+sy(2)*(IY-cy)+0.5,1);
      cnt = cnt+1;
    end
    %quadrant 2
    for w=1:nwpq
      FX{s}{cnt} = -ctr + (w-0.5)*2*ctr/nwpq;      FY{s}{cnt} = ctr;
      s2 = -FX{s}{cnt}/ctr;
      [nx,ny] = size(C{s}{cnt});
      NX{s}{cnt} = nx;      NY{s}{cnt} = ny;
      sx = [1/nx, s2*1/nx];      sy = [0, 1/ny];
      cx = nx/2+1;      cy = ny/2+1;
      [IX, IY] = ndgrid(1:nx, 1:ny);
      SX{s}{cnt} = 1 + m * mod(sx(1)*(IX-cx)+sy(1)*(IY-cy)+0.5,1);
      SY{s}{cnt} = 1 + n * mod(sx(2)*(IX-cx)+sy(2)*(IY-cy)+0.5,1);
      cnt = cnt+1;
    end
    %quadrant 3
    for w=1:nwpq
      FX{s}{cnt} = ctr;      FY{s}{cnt} = ctr - (w-0.5)*2*ctr/nwpq;
      s2 = FY{s}{cnt}/ctr;
      [nx,ny] = size(C{s}{cnt});
      NX{s}{cnt} = nx;      NY{s}{cnt} = ny;
      sx = [1/nx, 0];      sy = [-s2*1/ny, 1/ny];
      cx = nx/2+1;      cy = ny/2+1;
      [IX, IY] = ndgrid(1:nx, 1:ny);
      SX{s}{cnt} = 1 + m * mod(sx(1)*(IX-cx)+sy(1)*(IY-cy)+0.5,1);
      SY{s}{cnt} = 1 + n * mod(sx(2)*(IX-cx)+sy(2)*(IY-cy)+0.5,1);
      cnt = cnt+1;
    end
    %quadrant 4
    for w=1:nwpq
      FX{s}{cnt} = ctr - (w-0.5)*2*ctr/nwpq;      FY{s}{cnt} = -ctr;
      s2 = FX{s}{cnt}/ctr;
      [nx,ny] = size(C{s}{cnt});
      NX{s}{cnt} = nx;      NY{s}{cnt} = ny;
      sx = [1/nx, s2*1/nx];      sy = [0, 1/ny];
      cx = nx/2+1;      cy = ny/2+1;
      [IX, IY] = ndgrid(1:nx, 1:ny);
      SX{s}{cnt} = 1 + m * mod(sx(1)*(IX-cx)+sy(1)*(IY-cy)+0.5,1);
      SY{s}{cnt} = 1 + n * mod(sx(2)*(IX-cx)+sy(2)*(IY-cy)+0.5,1);
      cnt = cnt+1;
    end
    ctr = ctr*2;
  end
  
  %last
  SX{nbs} = cell(1,1);  SY{nbs} = cell(1,1);
  FX{nbs} = cell(1,1);  FY{nbs} = cell(1,1);
  NX{nbs} = cell(1,1);  NY{nbs} = cell(1,1);
  
  FX{nbs}{1} = 0;  FY{nbs}{1} = 0;
  [nx,ny] = size(C{nbs}{1});
  NX{nbs}{1} = nx;  NY{nbs}{1} = ny;
  sx = m/nx;  sy = n/ny;
  [TX,TY] = ndgrid(0:nx-1,0:ny-1);
  SX{nbs}{1} = 1 + TX*sx;  SY{nbs}{1} = 1 + TY*sy;
  
