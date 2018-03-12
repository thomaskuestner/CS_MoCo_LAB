function w = IteratedSineWindow(t)
% IteratedSine -- smooth window 
%  Usage
%    w = IteratedSineWindow(t)
%  Inputs
%    t     abscissa values for window evaluation
%  Outputs
%    w     beta(t + 1/2) beta(1/2 - t) 
%                beta is the iterated sine cut-off function
%                0 if t > 1 and 0 if t < -1.
%
% See Also IteratedSine, MakeWindow
%
% By Emmanuel Candes, 2003-2004

   
        w = IteratedSine(t+1/2).*IteratedSine(1/2 - t);
	
	
	
	