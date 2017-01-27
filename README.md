# CS_MoCo_LAB
Compressed Sensing and Motion Correction LAB: An MR acquisition and reconstruction system

Generate a Compressed Sensing (CS) accelerated MR sequence and reconstruct the acquired data online on the scanner by means of Gadgetron or offline on an external workstation.

## Acquistion
- Generic subsampling class for Compressed Sensing acquisitions
- CS accelerated gradient echo sequence (2D, 2D+time, 3D, 3D+time): CS_FLASH (Siemens, VB20P)
- MR motion imaging sequence (4D, 5D): CS_Retro (Siemens, VB20P) or refer to https://github.com/thomaskuestner/4DMRImaging

## Reconstruction
- CS and motion-resolved reconstruction system for Gadgetron (C++)
- CS reconstruction system in Matlab (including a GUI)

--------------------------------------------------------
Please read LICENSE file for licensing details.

Detailed information and installation instructions are available at:

https://sites.google.com/site/kspaceastronauts/compressed-sensing/cslab
