function [ y ] = operator_frame( x, n1,n2,n3, flag_complex, flag_in_trans, flag_out_trans, TransformOperator )
% represents base function for operator used
% Form Ax = y needed
% x; reshaped so A can be applied. (matrix and vector valid)
% y; reshaped into vector.
% flag_in_trans == 0: Transform takes standard input (mostly used)
% flag_in_trans == 1: Transform needs specific input (array)
% flag_out_trans == 0 : Transform provides standard output (depends on transform)
% flag_out_trans == 1: Transform forced to output vector (not in use)

% Transform: function handle, representing the operator

% initialize operator example:
% ST = @(x) sparse_trans(x,waveletStages,waveletFilterName_l1);
% ST_CG = @(x) operator_frame( x,n1,n2,flag_complex,flag_in_trans, ST )

% M. Fischer, April 2016

%%
if (flag_in_trans == 0)
    % x = reshape(x,1,n1*n2*n3); not needed atm
elseif (flag_in_trans == 1)
    if (isvector(x) && ~iscell(x))
        x = reshape(x,n1,n2,n3);
    end;
end;

if flag_complex
    if (~iscell(x))
        x_real = real(x);
        x_imag = imag(x);
    else
        x_real = cell(size(x,1),1);
        x_imag = x_real;
        for k = 1:size(x,1)
            x_real{k} = real(x{k});
            x_imag{k} = imag(x{k});
        end;
    end;

    y_real = TransformOperator(x_real);
    y_imag = TransformOperator(x_imag);
    
    if (~iscell(y_real))
        y = y_real + 1i*y_imag;
    else
        y = cell(size(y_real,1),1);
        for k = 1:size(y_real,1)
            y{k} = y_real{k} + 1i*y_imag{k};
        end;
    end;
else
    y = TransformOperator(x);
end;

if (flag_out_trans == 1)
    if (isvector(y) && ~iscell(x))
        y = reshape(y,n1,n2,n3);
    end;
end;

end

