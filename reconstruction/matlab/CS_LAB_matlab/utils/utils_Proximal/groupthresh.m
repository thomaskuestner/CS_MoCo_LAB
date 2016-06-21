function [x_thresh] = groupthresh(x,threshold,nCha,regularizer)
% thresholding for groups with same index (MMV case)
% atm l0, l1, p-shrinkage, MCP, SCAD and Arctangent and LOG penalty supported

% M. Fischer, April 2016

%% prelimanary test values
a = 0.5/(nCha*threshold); % a = 1 was way better
p = -0.5; % p-shrinkage; often used with {-0.5,0,0.5}
x_thresh = cell(1,nCha);
%% calculate l2-norm of groups:
x_group = zeros(size(x{1,1}));
for j=1:nCha
    x_group = x_group + abs(x{1,j}).^2;
end;
x_group = sqrt(x_group);

%% thresholding:
switch regularizer
    case 'l0'
        x_group_thresh = wthresh(x_group,'h',threshold);
        % x_group_thresh = 1./(x_group - threshold) .* max(x_group - threshold, 0);
    case 'lp'
        warning('Regularizer has not been implemented');
        % u may take a look at GISA.
    case 'l1'
        x_group_thresh = softthresh_real(x_group,threshold);
        %x_group_thresh = wthresh(x_group,'s',threshold);
    case 'MCP'
        for l = 1:size(x_group,5)
            for k = 1:size(x_group,4)
                for j = 1:size(x_group,3)
                    if size(x_group,5) == 1 && size(x_group,4) == 1 && size(x_group,3) == 1 % case of 1xN vector
                        x_group_thresh(:,1,j,k,l) = proximalRegC(permute(x_group,[2 1 3 4 5]), size(x_group(1,:,j,k,l),2), threshold, 2, 4); % from GIST, theta >= 2
                        x_group_thresh = permute(x_group_thresh,[2 1 3 4 5]);
                    else
                        for h = 1:size(x_group,2) % multidimensional case
                            x_group_thresh(:,h,j,k,l) = proximalRegC(x_group(:,h,j,k,l), size(x_group(:,h,j,k,l),1), threshold, 2, 4); % from GIST, theta >= 2
                            % x_group_thresh = x_group_thresh';
                        end;
                    end;
                end;
            end;
        end;
    case 'SCAD'
        for l = 1:size(x_group,5)
            for k = 1:size(x_group,4)
                for j = 1:size(x_group,3)                       
                    if size(x_group,5) == 1 && size(x_group,4) == 1 && size(x_group,3) == 1 % case of 1xN vector
                        x_group_thresh(:,1,j,k,l) = proximalRegC(permute(x_group,[2 1 3 4 5]), size(x_group(1,:,j,k,l),2), threshold, 3, 3); % from GIST, theta >= 3
                        x_group_thresh = permute(x_group_thresh,[2 1 3 4 5]);
                    else
                        for h = 1:size(x_group,2) % multidimensional case
                            x_group_thresh(:,h,j,k,l) = proximalRegC(x_group(:,h,j,k,l), size(x_group(:,h,j,k,l),1), threshold, 3, 3); % from GIST, theta >= 3
                            % x_group_thresh = x_group_thresh';
                        end;
                    end;
                end;
            end;
        end;
    case 'log'
        x_group_thresh = (x_group/2-1/(2*a) + sqrt(x_group./2 + 1/(2*a)).^2 - threshold/a); % from ParekhDenoising
    case 'ATAN'
        x_group_thresh = atanT2(x_group,threshold,a); % from ParekhDenoising
    case 'PSHRINK'
        x_group_thresh = max(x_group - threshold^(2-p)* x_group.^(p-1),0); % from "IPS" paper
    case 'FIRM'
        warning('Regularizer has not been implemented');
        % use MCP instead.
    otherwise
        warning('No valid regularizer chosen');
end;

for j=1:nCha
    x_thresh{1,j} = x{1,j} ./x_group .* x_group_thresh;
end;

