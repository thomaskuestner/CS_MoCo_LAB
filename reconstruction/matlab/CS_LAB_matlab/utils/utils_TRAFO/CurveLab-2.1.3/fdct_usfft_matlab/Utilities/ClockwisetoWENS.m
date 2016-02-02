function D = ClockwisetoWENS(Cs);
% WENStoClockwise: Converts a data structure where wedge
%                  start from top left wedge and increasing in a
%                  clockwise fashion, into another structure with
%                  orientations organized by cardinal points. 
%
% See Also
%   SeparateAngles, Adj_DetailCurveCoeff, WENStoClockwise
     
  
    nw = length(Cs)/4;    [nl,ns] = size(Cs{1});
    D = zeros(4,nw,nl,ns);
    cnt = 1;
    for w=1:nw
      D(3,w,:,:) = Cs{cnt}; cnt = cnt+1;
    end
    for w=nw:-1:1
      D(2,w,:,:) = Cs{cnt}.'; cnt = cnt+1;
    end
    for w=1:nw
      D(4,w,:,:) = Cs{cnt}; cnt = cnt+1;
    end
    for w=nw:-1:1
      D(1,w,:,:) = Cs{cnt}.'; cnt = cnt+1;
    end

   


  