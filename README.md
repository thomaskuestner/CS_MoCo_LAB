# CS_LAB
Compressed Sensing LAB: An MR acquisition and reconstruction system

Generate a Compressed Sensing (CS) accelerated MR sequence and reconstruct the acquired data online on the scanner by means of Gadgetron or offline on an external workstation.

## Acquistion
- Generic subsampling class for Compressed Sensing acquisitions
- precompiled CS accelerated gradient echo sequence (2D, 2D+time, 3D, 3D+time): CS_FLASH (Siemens, VB20P)
- for a full 4D (3D+time) CS accelerated gradient echo sequence under free movement, please refer to: CS_Retro in https://github.com/thomaskuestner/4DMRImaging

## Reconstruction
- CS reconstruction system for Gadgetron (C++)
- CS reconstruction system in Matlab (including a GUI)

--------------------------------------------------------
Please read LICENSE file for licensing details.

Detailed information and installation instructions are available at:

https://sites.google.com/site/kspaceastronauts/compressed-sensing/cslab
