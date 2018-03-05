function [ y ] = trafo_WaT( x,wavelevel,wavename,G_mat,Gt_mat )
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

    [x_e1, S] = wavedec2(x,wavelevel,wavename);
    x_e2 = G_mat*x_e1'; 
    y = waverec2(x_e1+(Gt_mat*x_e2)',S,wavename);

end

