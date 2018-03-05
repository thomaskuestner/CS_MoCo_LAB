function xrec=DoPyrRec(L)
Pyr_mode=1.5;
smooth_func = @rcos;
xrec = PyrNDRec_mm(L, 'S', Pyr_mode, smooth_func);