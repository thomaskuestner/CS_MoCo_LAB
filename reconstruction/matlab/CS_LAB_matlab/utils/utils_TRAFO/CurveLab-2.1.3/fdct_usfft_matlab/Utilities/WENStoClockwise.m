function Cs = WENStoClockwise(D);
% WENStoClockwise: Converts a data structure where wedge orientations
%                  are organized by cardinal point (in the order West
%                  East North South) into another structure with
%                  orientations starting from top left wedge and
%                  increasing in a clockwise fashion. 
%
% See Also
%   SeparateAngles, DetailCurveCoeff, ClockwisetoWENS
  

   [nq,nw,nl,ns] = size(D);
   Cs = cell(1,nq*nw);
   cnt = 1;

    for w = 1:nw
      Cs{cnt} = squeeze(D(3,w,:,:)); cnt = cnt+1;
    end
    for w=nw:-1:1
      Cs{cnt} = squeeze(D(2,w,:,:)).'; cnt = cnt + 1;
    end
    for w=1:nw
      Cs{cnt} = squeeze(D(4,w,:,:)); cnt = cnt+1;
    end
    for w=nw:-1:1
      Cs{cnt} = squeeze(D(1,w,:,:)).'; cnt = cnt + 1;
    end
   


  