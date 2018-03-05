--- prerequisites ---
- Armadillo 6
- blas resp. openblas (dev)
- lapack (dev)
- Insight Segmentation and Registration Toolkit (ITK)

- for parallel computing: compiler with openmp standard (recommended)

--- (standalone) installation instruction ---
mkdir /home/LAP
cd /home/LAP
cmake /PATH_TO_SRC
make

--- testing ---
call programm with image paths (first image, second image, result image) and image sizes (rows, columns, slices)