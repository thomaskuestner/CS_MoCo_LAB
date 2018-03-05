function [ y ] = MTV_2D( x,lambda,n1,n2 )
% solves the anisotrope TV-proximal operator
% using the 1D-TV Code by Condat
%
% (c) Marc Fischer, Thomas Kuestner
% ---------------------------------------------------------------------

if (0 < lambda) && (lambda < 1) % temporary check

    yt_x = zeros(n1,n2);
    yt_y = yt_x;

    vec_y = cell(n1,1);
    vec_x = cell(n2,1);

    % uses PSA (just uses the mean)
    % no input validation (see C Source Mex-File)
    
            for i=1:n2
                vec_x{i,1} = x(:,i)';
                vec_x{i,1} = MTV_1D(double(lambda),double(vec_x{i,1}));
                yt_x(:,i) = vec_x{i,1}';
            end;
            for l=1:n1
                vec_y{l,1}= x(l,:);
                vec_y{l,1} = MTV_1D(double(lambda),double(vec_y{l,1}));
                yt_y(l,:) = vec_y{l,1};
            end;

            y = (yt_x + yt_y) ./2;
            y = single(y);
          
else
    y = zeros(n1,n2);
end

