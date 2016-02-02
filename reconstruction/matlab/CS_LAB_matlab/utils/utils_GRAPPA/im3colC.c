#include <mex.h>
#include <matrix.h>
// rearrange column vector to 3D array (complex-valued)
//
// (c) Thomas Kuestner
// ---------------------------------------------------------------------

// stupid Matlab mex is unable to compile with if statements checking boolean variables
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  const mxArray *I = prhs[0];
  double *indataReal = mxGetPr(I);
  //double *indataImag;
  //if( mxIsComplex(I) )
  double *indataImag = mxGetPi(I);
  double *patchSize = mxGetPr(prhs[1]);
  const int *size = mxGetDimensions(I);
  int J = (int)patchSize[0];
  int K = (int)patchSize[1];
  int H = (int)patchSize[2];
  int M = size[0];
  int N = size[1];
  int P = size[2];

  int numPatches = (M - J + 1)*(N - K + 1)*(P - H + 1);
  int out_rows = J*K*H;
  int out_cols = numPatches;
  //const mxArray *isComplex;
  //double *outdataReal;
  //double *outdataImag;
  mxArray *out = mxCreateDoubleMatrix( out_rows, out_cols, mxCOMPLEX );
  double *outdataReal = mxGetPr(out);
  double *outdataImag = mxGetPi(out);
  /*if ( !mxIsComplex(I) ) {
    out = mxCreateDoubleMatrix( out_rows, out_cols, mxREAL );  
    outdataReal = mxGetPr(out);
    isComplex = mxCreateLogicalScalar(false);
  } else {
    out = mxCreateDoubleMatrix( out_rows, out_cols, mxCOMPLEX );
    outdataImag = mxGetPi(out);
    isComplex = mxCreateLogicalScalar(true); 
  }*/
  
  
  int patch = 0;
  int h_offset = 0;
  for( h_offset = 0; h_offset < P-H+1; h_offset++ ){
	int k_offset = 0;
    for( k_offset = 0; k_offset < N-K+1; k_offset++ ){
	  int j_offset = 0;
      for( j_offset = 0; j_offset < M-J+1; j_offset++ ){
        int row = 0;
		int h = 0;
        for( h = 0; h < H; h++ ){
		  int k = 0;
          for( k = 0; k < K; k++ ){
		    int j = 0;
            for( j = 0; j < J; j++ ){
                //if(mxIsLogicalScalarTrue(isComplex)) {
                    outdataReal[patch*out_rows + row] = 
                indataReal[ (j_offset+j) + (k_offset+k)*M + (h_offset+h)*M*N ];
                    outdataImag[patch*out_rows + row] = 
                indataImag[ (j_offset+j) + (k_offset+k)*M + (h_offset+h)*M*N ];
                /*} else {
                    outdataReal[patch*out_rows + row] = 
                indataReal[ (j_offset+j) + (k_offset+k)*M + (h_offset+h)*M*N ];
                }*/
              ++row;
            }}}
      ++patch;
      }}}
  plhs[0] = out;
}