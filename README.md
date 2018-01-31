# CS_MoCo_LAB [![Build Status](https://semaphoreci.com/api/v1/thomaskuestner/cs_moco_lab/branches/master/shields_badge.svg)](https://semaphoreci.com/thomaskuestner/cs_moco_lab)
Compressed Sensing and Motion Correction LAB: An MR acquisition and reconstruction system

Generate a Compressed Sensing (CS) accelerated MR sequence and reconstruct the acquired data online on the scanner by means of Gadgetron or offline on an external workstation.

## Acquistion
- Generic subsampling class for Compressed Sensing acquisitions
- CS accelerated gradient echo sequence (2D, 2D+time, 3D, 3D+time): CS_FLASH (Siemens, VB20P)
- MR motion imaging sequence (4D, 5D): CS_Retro (Siemens, VB20P, VE11C, VE11P)

## Reconstruction
- CS and motion-resolved reconstruction system for Gadgetron (C++)
- CS reconstruction system in Matlab (including a GUI)

## Applications
### Compressed Sensing MRI
https://sites.google.com/site/kspaceastronauts/compressed-sensing/espresso <br/>
https://sites.google.com/site/kspaceastronauts/compressed-sensing/cslab
### Motion MR imaging
https://sites.google.com/site/kspaceastronauts/motion-correction/4d-mr-imaging <br/>
https://github.com/thomaskuestner/4DMRImaging
### PET/MR motion correction
https://sites.google.com/site/kspaceastronauts/motion-correction/pet-mr-motion-correction
- acquisition: CS_Retro (Siemens, VB20P, VE11C, VE11P)
- reconstruction: CS_LAB in gadgetron + data emitters and injectors
### Image Registration
https://sites.google.com/site/kspaceastronauts/motion-correction/mocogui

## Installation
More information: https://sites.google.com/site/kspaceastronauts/compressed-sensing/cslab#TOC-Prerequisites
### Gadgetron ![Gadetron Status1](https://img.shields.io/badge/Gadgetron-v2.5.0-brightgreen.svg) ![Gadgetron Status2](https://img.shields.io/badge/Gadgetron-v3.15.0-red.svg)
**standalone**<br/>
```
git clone https://github.com/thomaskuestner/CS_MoCo_LAB.git
mkdir build
cd build
cmake ../
make
sudo make install
git clone https://github.com/thomaskuestner/CS_MoCo_LAB.git
cd CS_MoCo_LAB/gadgetron/CS_LAB_Gadget
make
```

**docker**
```
docker pull kspaceastronauts/cs_moco_lab
docker run -it --volume $(pwd):/opt/data kspaceastronauts/cs_moco_lab
```

### Matlab
```
git clone https://github.com/thomaskuestner/CS_MoCo_LAB.git
cd CS_MoCo_LAB/matlab/CS_LAB_GUI
matlab CS_LAB_GUI
```

--------------------------------------------------------
Please read LICENSE file for licensing details.

Detailed information and installation instructions are available at: <br/>
https://sites.google.com/site/kspaceastronauts/compressed-sensing/cslab
