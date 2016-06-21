//////////////////////////////////////////////////////////////////////////
//
//	Pooran Singh Negi and Demetrio Labate
//
//	Department of Mathematics
//	University of Houston
//
//////////////////////////////////////////////////////////////////////////
//
//	This software is provided "as-is", without any express or implied
//	warranty. In no event will the authors be held liable for any 
//	damages arising from the use of this software.
//
//////////////////////////////////////////////////////////////////////////

Version 1.0 4/30/13

---------------------------------------------------------------------------------------
INSTALL INSTRUCTIONS


For running 3DShearlet code user need to  run
convnfft_install.m from CONVNFFT_Folder  first.
---------------------------------------------------------------------------------------


---------------------------------------------------------------------------------------
Restriction on the size of .mat file

Due to nature of upsampling and downsampling for better denoising performance
dimension of data need to be divisible by 3*2^(number of decomposition level required).
---------------------------------------------------------------------------------------



---------------------------------------------------------------------------------------
RUNNING DENOISING DEMO

Gaussian Denoising script is  in DenoiseDemo folder with file name main.m
If user need to save data please use aproriate save statement after deonising is done.
Currrent script setting  denoises coastguard_sequence with simulated noise.
This script can be run in background in Unix/Linux using hmatbg in DenoiseDemo folder.
./hmatbg main.m outfile 
hmatbg can be modified  to change priority level using nice command.


For data corrupted by unknown Gaussian noise parameter, sdest3 funtion in Util directory
can used to estimate standard deviation of the noise which is based on median of the wavelet coefficients at the finer scale.
---------------------------------------------------------------------------------------

---------------------------------------------------------------------------------------
NOTE: Files in 3DBP directory for doing bandpass are from SurfBox toolbox implemented by  Yue Lu and Minh N. Do
---------------------------------------------------------------------------------------
