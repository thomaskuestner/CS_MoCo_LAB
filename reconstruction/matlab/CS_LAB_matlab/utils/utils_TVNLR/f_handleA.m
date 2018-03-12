function y = f_handleA(A,u,mode,opts)
num_rows = opts.row;
num_cols = opts.col;
block_size = opts.block_size;
ratio = opts.ratio;
m = round(ratio*block_size*block_size);
switch mode
    case 1
        %y = A*u;
        x = reshape(u,num_rows,num_cols);
        x_hat = im2col(x, [block_size block_size], 'distinct');
        y = A*x_hat;
        y = y(:);
    case 2
        %y = (u'*A)';
        y = reshape(u,m,length(u)/m);
        x_hat = A'*y;
        x = col2im(x_hat, [block_size block_size], ...
            [num_rows num_cols], 'distinct');
        y = reshape(x,num_rows*num_cols,1);
    otherwise
        error('Unknown mode passed to f_handleA in ftv_cs.m');
end

end