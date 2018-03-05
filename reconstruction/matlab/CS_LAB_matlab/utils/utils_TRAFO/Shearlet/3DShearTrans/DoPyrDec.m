function L= DoPyrDec(x,level)

%% use the raised-cosine function to get a smooth transition band. 
smooth_func = @rcos;
Pyr_mode=1.5;
L = PyrNDDec_mm(x, 'S', level, Pyr_mode, smooth_func);