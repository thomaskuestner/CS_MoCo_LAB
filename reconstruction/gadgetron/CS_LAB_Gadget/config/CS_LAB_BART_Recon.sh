#!/bin/bash

# List of available default parameters
# recon_matrix_x;
# recon_matrix_y;
# recon_matrix_z;
# FOV_x;
# FOV_y;
# FOV_z;
# acc_factor_PE1;
# acc_factor_PE2;
# reference_lines_PE1;
# reference_lines_PE2;


debug=false

if "$debug"; then
	set -x
fi

bart ecalib -S input_data maps
bart pics -R Q:.01 -S input_data maps img_recon0
bart slice 4 0 img_recon0 ims

if "$debug";then
	set +x
fi
