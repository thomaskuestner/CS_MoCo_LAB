function [ y ] = sparse_trans( x, TransformSpecifics, reconDIM )
% used to set up transform handle

% M. Fischer, April 2016

%%

if ~iscell(x) && isreal(x) || iscell(x) && isreal(x{1}) % is this error safe?
    if strcmp(reconDIM,'2D')
        switch TransformSpecifics.TransformName
            case 'dwt'
                y = wavedec2(x,TransformSpecifics.Stages,TransformSpecifics.FilterName);
            case 'cplxdt_matlab' % matlab variant:
                y_helper = dddtree2('cplxdt',x,3,TransformSpecifics.dfFilter{1},TransformSpecifics.dfFilter{2});
                y = y_helper.cfs;
%                 y = cell(TransformSpecifics.Stages+1,1);
%                 for k = 1:TransformSpecifics.Stages+1 % change to array to apply thresholding
%                     y{k} = y_helper.cfs{k}; % cfs{i} contains a 5D array.
%                 end;
            case 'cplxdt' % NOTE! n1 & n2 must be divisible by 2^J.
                if TransformSpecifics.flagExtend
                    x = img_extend(x,TransformSpecifics.extend_x,TransformSpecifics.extend_y,0);
                end;
                y_helper = cplxdual2D(x,TransformSpecifics.Stages,TransformSpecifics.Faf,TransformSpecifics.af);
                y = cell(TransformSpecifics.Stages+1,1);
                for l = 1:TransformSpecifics.Stages % wavletStages (j)
                    for k = 1:2 % real & imag (i)
                        for j = 1:2 % orientation#1 (d1)
                            for h = 1:3 % orientation#2 (d2)
                                y{l}(:,:,h,j,k) = y_helper{1,l}{1,k}{1,j}{1,h}; % contains a 5D array.
                            end;
                        end;
                    end;
                end;
                for j = 1:2 % orientation#1 (d1)
                    for h = 1:2 % orientation#2 (d2)
                        y{TransformSpecifics.Stages+1}(:,:,h,j,1) = y_helper{1,TransformSpecifics.Stages+1}{1,j}{1,h};
                    end;
                end;
            case 'shearlet'
                if TransformSpecifics.flagExtend
                    x = img_extend(x,TransformSpecifics.extend_x,TransformSpecifics.extend_y,0);
                end;
                [y,~]= nsst_dec2(x,TransformSpecifics.shear_parameters,TransformSpecifics.lpfilt);
                y = y.';
            case 'shearlet_nontight'
                y = SLsheardec2D(x,TransformSpecifics.shearletSystem);
            otherwise
                warning('No correct Transform chosen');
        end;
    elseif strcmp(reconDIM,'3D')
        switch TransformSpecifics.TransformName
            case 'dwt'
                y_helper = wavedec3(x,TransformSpecifics.Stages,TransformSpecifics.FilterName);
                y = y_helper.dec;
            case 'cplxdt_matlab' % matlab variant:
                warning('cplxdt_matlab does not exist');
            case 'cplxdt' % NOTE! n1 & n2 must be divisible by 2^J.
                if TransformSpecifics.flagExtend
                    x = img_extend(x,TransformSpecifics.extend_x,TransformSpecifics.extend_y,TransformSpecifics.extend_z);
                end;
                y_helper = cplxdual3D(x,TransformSpecifics.Stages,TransformSpecifics.Faf,TransformSpecifics.af);
                y = cell(TransformSpecifics.Stages+1,1);
                for l = 1:TransformSpecifics.Stages % wavletStages (j)
                    for k = 1:2 % (m)
                        for j = 1:2 % (n)
                            for h = 1:2 % (p)
                                for g = 1:7 % (d)
                                    y{l}(:,:,:,g,h,j,k) = y_helper{1,l}{1,k}{1,j}{1,h}{1,g}; % contains a 7D array.
                                end;
                            end;
                        end;
                    end;
                end;
                for k = 1:2 % (m)
                    for j = 1:2 % (n)
                        for g = 1:2 % (d) not 1:7
                        y{TransformSpecifics.Stages+1}(:,:,:,g,1,j,k) = y_helper{1,TransformSpecifics.Stages+1}{1,k}{1,j}{1,g};
                        end;
                    end;
                end;
            case 'shearlet'
                if TransformSpecifics.flagExtend
                    x = img_extend(x,TransformSpecifics.extend_x,TransformSpecifics.extend_y,TransformSpecifics.extend_z);
                end;
                y_helper = ShearTransform3D(x,TransformSpecifics.shearingFilter);
                decomp_cells = 3*(4*4+2*2); % (fixed in get_transformSpecifics) dims:3, level:2 ,directions: 4x4 & 2x2 
                % hard coded solution: (feel free to replace it)
                y = cell(1+decomp_cells,1);
                i = 1;
                for l = 1:3 % dims
                    k = 1; % (level)
                    for j = 1:4 % (direction)
                        for h = 1:4 % (direction)
                            y{i} = y_helper.D{l,k}{j,h};
                            i = i+1;
                        end;
                    end;
                    k =2;
                    for j = 1:2 % (direction)
                        for h = 1:2 % (direction)
                            y{i} = y_helper.D{l,k}{j,h};
                            i = i+1;
                        end;
                    end;
                end;
                y{i+1} = y_helper.A;
            case 'shearlet_nontight'
                if TransformSpecifics.flagExtend
                    x = img_extend(x,0,0,TransformSpecifics.extend_z);
                end;
                y = SLsheardec3D(x,TransformSpecifics.shearletSystem);
            otherwise
                warning('No correct Transform chosen');
        end;
    else
        error('ReconDim is not supported');
    end;   
else % should always be real if operator_frame.m is used.
    error('use operator_frame.m');
end;

end

