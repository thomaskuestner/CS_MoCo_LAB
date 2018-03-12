function w = MakeSineWindow(boxlen,sigma)
% MakeSineWindow -- Make a nice window out of the IteratedSine function.
%  Inputs
%    boxlen   abscissa values for window evaluation -boxlen <= t < boxlen 
%    sigma    scale of the window, controls the size of the
%             effective support           
%  Outputs
%    w        values of the window
%  See Also
%    IteratedSine
%
% By Emmanuel Candes, 2003-2004
  
  ix = ((-boxlen:(boxlen -1)) + .5)./boxlen;
  w = IteratedSineWindow(ix./sigma);
