%% Export checkerboard filters to files
%% Used by Surfbox-C++

%% First created: 04-09-2006
%% Last modified: 04-09-2006

function export_cbd_filters(bo)

fname = ['cbd_coeffs_bo_' num2str(bo) '.surf'];
fid = fopen(fname, 'w');

for n = 1 : 2
    
    if (n == 1)
        %% decomposition filters
        flt = get_cbd_filters_load(bo, 'd');
    else
        %% reconstruction filters
        flt = get_cbd_filters_load(bo, 'r');
    end

    %% filter 1
    sz = fliplr(size(flt{1}));
    fwrite(fid, sz, 'int');
    center = fliplr(flt{2}) - 1;
    fwrite(fid, center, 'int');

    coeffs = flt{1};
    fwrite(fid, coeffs, 'double');

    %% filter 2
    sz = fliplr(size(flt{3}));
    fwrite(fid, sz, 'int');
    center = fliplr(flt{4}) - 1;
    fwrite(fid, center, 'int');

    coeffs = flt{3};
    fwrite(fid, coeffs, 'double');
end

fclose(fid);