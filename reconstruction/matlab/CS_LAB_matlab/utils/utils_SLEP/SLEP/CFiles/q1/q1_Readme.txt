This folder contains the following C files:

epph.h is the head file that contains the most important subroutines.

------------------------------------------------------------------------
ep1R is for solving the problem 
   min  1/2 ( \|x- u\|^2 +  \|t-v\|^2 )
   s.t. |x_j| <= t_j
  
   x, u, t, and v are of size nx1
   x_j and t_j denotes the j-th element of x and t, respectively
------------------------------------------------------------------------
   
------------------------------------------------------------------------
ep21R is for solving the problem 
   min  1/2 ( \|x- u\|^2 +  \|t-v\|^2 )
   s.t. \|x_j\|_2 <= t_j

   x and u are of size n x k
   t and v are of size nx1
   x_j denotes the j-th row of x
   t_j denotes the j-th element of t
------------------------------------------------------------------------

------------------------------------------------------------------------
ep21d is for solving the problem 
   min  1/2 \|x- v\|^2 
   s.t. \sum_j \|x_j\|_2 <= z

   x and v are of size n x k
   x_j denotes the j-th row of x
------------------------------------------------------------------------

------------------------------------------------------------------------
eplb is for solving the problem 
   min  1/2 \|x- v\|^2 
   s.t. \|x\|_1 <= z

   x and v are of size n x 1
   \|x\|_1 denotes 1-norm of x
------------------------------------------------------------------------

------------------------------------------------------------------------
eppMatrix is for solving the problem 
   min  1/2 ||x- v||_2^2 + rho * \sum_j ||x_j||_p 

   x and v are of size n x k
   x_j denotes the j-th row of x
------------------------------------------------------------------------

------------------------------------------------------------------------
eppVector is for solving the problem 
   min  1/2 ||x- v||_2^2 + rho * \sum_j gWeight_j * ||x_j||_p 

   x and v are of size n x 1
   gWeight is of size k x 1
   x_j denotes the j-th group of x (the entries are adjacent)
   gWeight_j is the weight for the j-th group
------------------------------------------------------------------------

------------------------------------------------------------------------
eppVectorR is for solving the problem 
   min  1/2 ( ||x- u||_2^2 + ||t-v||_2^2 )
   s.t.  ||x_j||_2 <= t_j

   x and u are of size n x 1
   t and v are of size k x 1
   x_j denotes the j-th group of x (the entries are adjacent)
------------------------------------------------------------------------