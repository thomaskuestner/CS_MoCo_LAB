function [ y ] = back_trans( x, TransformSpecifics, reconDIM, n1, n2, n3 )
% used to set up transform handle

% M. Fischer, April 2016

%%

if ~iscell(x) && isreal(x) || iscell(x) && isreal(x{1}) % is this error safe?
    if strcmp(reconDIM,'2D')
        switch TransformSpecifics.TransformName
            case 'dwt'
                y = waverec2(x,TransformSpecifics.waveS,TransformSpecifics.FilterName);
            case 'cplxdt_matlab' % matlab variant:
                x_helper = TransformSpecifics.wt_struct;
%                 for k = 1:TransformSpecifics.Stages+1 % change to array to apply thresholding
%                     x_helper.cfs{k} = x{k}; % cfs{i} contains a 5D array.
%                 end;
                x_helper.cfs = x;
                y = idddtree2(x_helper);
            case 'cplxdt' % NOTE! n1 & n2 must be divisible by 2^J.
                for l = 1:TransformSpecifics.Stages % wavletStages (j)
                    for k = 1:2 % real & imag (i)
                        for j = 1:2 % orientation#1 (d1)
                            for h = 1:3 % orientation#2 (d2)
                                x_helper{1,l}{1,k}{1,j}{1,h} = x{l}(:,:,h,j,k); % contains a 5D array.
                            end;
                        end;
                    end;
                end;
                for j = 1:2 % orientation#1 (d1)
                    for h = 1:2 % orientation#2 (d2)
                        x_helper{1,TransformSpecifics.Stages+1}{1,j}{1,h} = x{TransformSpecifics.Stages+1}(:,:,h,j,1);
                    end;
                end;
                y = icplxdual2D(x_helper, TransformSpecifics.Stages, TransformSpecifics.Fsf, TransformSpecifics.sf);
                y = y(1:n1,1:n2);
            case 'shearlet'
                y = nsst_rec2(x.',TransformSpecifics.shear_f,TransformSpecifics.lpfilt);
                y = y(1:n1,1:n2);
            case 'shearlet_nontight'
                y = SLshearrec2D(x,TransformSpecifics.shearletSystem);
            otherwise
                warning('No correct Transform chosen');
        end;    
    elseif strcmp(reconDIM,'3D')
        switch TransformSpecifics.TransformName
            case 'dwt'
                x_helper = TransformSpecifics.WDEC;
                x_helper.dec = x;
                y = waverec3(x_helper);
            case 'cplxdt'
                for l = 1:TransformSpecifics.Stages % wavletStages (j)
                    for k = 1:2 % (m)
                        for j = 1:2 % (n)
                            for h = 1:2 % (p)
                                for g = 1:7 % (d)
                                    x_helper{1,l}{1,k}{1,j}{1,h}{1,g} = x{l}(:,:,:,g,h,j,k); % contains a 7D array.
                                end;
                            end;
                        end;
                    end;
                end;
                for k = 1:2 % (m)
                    for j = 1:2 % (n)
                        for g = 1:2 % (d) not 1:7
                        x_helper{1,TransformSpecifics.Stages+1}{1,k}{1,j}{1,g} = x{TransformSpecifics.Stages+1}(:,:,:,g,1,j,k);
                        end;
                    end;
                end;
                y = icplxdual3D(x_helper, TransformSpecifics.Stages, TransformSpecifics.Fsf, TransformSpecifics.sf);
                y = y(1:n1,1:n2,1:n3);
            case 'shearlet'
                decomp_cells = 3*(4*4+2*2); % (fixed in get_transformSpecifics) dims:3, level:2 ,directions: 4x4 & 2x2 
                % hard coded solution: (feel free to replace it)
                y = cell(1+decomp_cells,1);
                i = 1;
                for l = 1:3 % dims
                    k = 1; % (level)
                    for j = 1:4 % (direction)
                        for h = 1:4 % (direction)
                            x_helper.D{l,k}{j,h} = x{i};
                            i = i+1;
                        end;
                    end;
                    k =2;
                    for j = 1:2 % (direction)
                        for h = 1:2 % (direction)
                            x_helper.D{l,k}{j,h} = x{i};
                            i = i+1;
                        end;
                    end;
                end;
                x_helper.A = x{i+1};
                y = InvShearTransform3D(x_helper);
            case 'shearlet_nontight'
                if TransformSpecifics.flagExtend
                    x = img_extend(x,0,0,TransformSpecifics.extend_z);
                end;
                y = SLshearrec3D(x,TransformSpecifics.shearletSystem);
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