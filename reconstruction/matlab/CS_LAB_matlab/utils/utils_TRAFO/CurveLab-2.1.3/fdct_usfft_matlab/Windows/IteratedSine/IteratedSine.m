function w = IteratedSine(t)
% IteratedSine -- window for iterated sine cut-off function 
%  Usage
%    w = IteratedSine(t)
%  Inputs
%    t     abscissa values for window evaluation
%  Outputs
%    w     sin(pi/4(1 + sin(pi t))) 
%           1 if t > 1/2 and 0 if t < -1/2.
%  See Also
%    IteratedSineWindow, MakeSineWindow
%
% By Emmanuel Candes, 2003-2004


        phase = 1 + sin(pi*t);
	w     = sin(pi/4 * phase);
	
	t0 = find(t <= -1/2);
	if length(t0) > 0,
		w(t0) = zeros(1,length(t0));
	end
	t1 = find(t >= 1/2);
	if length(t1) > 0,
		w(t1) = ones(1,length(t1));
	end

	