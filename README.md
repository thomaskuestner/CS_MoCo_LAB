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
### Gadgetron ![Gadetron Status1](https://img.shields.io/badge/Gadgetron-v2.5.0-brightgreen.svg) ![Gadgetron Status2](https://img.shields.io/badge/Gadgetron-v3.15.0-green.svg)
**standalone**<br/>
```
git clone https://github.com/thomaskuestner/CS_MoCo_LAB.git
mkdir build
cd build
cmake ../
make
sudo make install
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

## References
<b>MR-based respiratory and cardiac motion correction for PET imaging</b>.<br/>
<i>Medical Image Analysis</i>, 2017.
<br/>Thomas K&uuml;stner, Martin Schwartz, Petros Martirosian, Sergios Gatidis, Ferdinand Seith, Christopher Gilliam, Thierry Blu, Hadi Fayad, Dimitris Visvikis, F. Schick, B. Yang, H. Schmidt and N.F Schwenzer.
<br/><a href="http://dx.doi.org/10.1016/j.media.2017.08.002">[doi]</a>&nbsp;<a onclick="toggleBibtex('thomaskuestner', '373a2d347e30bfae4639dd40daaba4f3', 'https://puma.ub.uni-stuttgart.de/bibtex/2373a2d347e30bfae4639dd40daaba4f3/thomaskuestner?format=bibtex'); return false;" href="https://puma.ub.uni-stuttgart.de/bibtex/2373a2d347e30bfae4639dd40daaba4f3/thomaskuestner?format=bibtex">[BibTeX]</a>&nbsp;<a onclick="toggleEndnote('thomaskuestner', '373a2d347e30bfae4639dd40daaba4f3', 'https://puma.ub.uni-stuttgart.de/bibtex/2373a2d347e30bfae4639dd40daaba4f3/thomaskuestner?format=endnote'); return false;" href="https://puma.ub.uni-stuttgart.de/bibtex/2373a2d347e30bfae4639dd40daaba4f3/thomaskuestner?format=endnote">[Endnote]</a>


<b>Compressed Sensing LAB: An MR acquisition and reconstruction system</b>. <br/>
<i>Proceedings of the ISMRM Workshop on Data Sampling and Reconstruction</i>. Sedona, AZ, USA, 2016.
<br/>Thomas K&uuml;stner, Martin Schwartz, Christian W&uuml;rslin, Petros Martirosian, Nina F. Schwenzer, Bin Yang and Holger Schmidt.
<br/><a onclick="toggleBibtex('thomaskuestner', 'a2f614af9059696b743e4a8cd0b1f7e4', 'https://puma.ub.uni-stuttgart.de/bibtex/2a2f614af9059696b743e4a8cd0b1f7e4/thomaskuestner?format=bibtex'); return false;" href="https://puma.ub.uni-stuttgart.de/bibtex/2a2f614af9059696b743e4a8cd0b1f7e4/thomaskuestner?format=bibtex">[BibTeX]</a>&nbsp;<a onclick="toggleEndnote('thomaskuestner', 'a2f614af9059696b743e4a8cd0b1f7e4', 'https://puma.ub.uni-stuttgart.de/bibtex/2a2f614af9059696b743e4a8cd0b1f7e4/thomaskuestner?format=endnote'); return false;" href="https://puma.ub.uni-stuttgart.de/bibtex/2a2f614af9059696b743e4a8cd0b1f7e4/thomaskuestner?format=endnote">[Endnote]</a>&nbsp;<div id="abs_a2f614af9059696b743e4a8cd0b1f7e4thomaskuestner" style="display:none;border:1px dotted grey;">

<b>Image Reconstruction System for Compressed Sensing Retrospective Motion Correction for the Application in Clinical Practice</b>. <br/>
<i>Proceedings of the International Society for Magnetic Resonance in Medicine (ISMRM)</i>. Singapore, 2016.
<br/>Martin Schwartz, Thomas K&uuml;stner, Christian W&uuml;rslin, Petros Martirosian, Nina F. Schwenzer, Fritz Schick, Bin Yang and Holger Schmidt.
<br/><a onclick="toggleBibtex('thomaskuestner', 'd4738bfd3b3a4eeefd763dec15715d9d', 'https://puma.ub.uni-stuttgart.de/bibtex/2d4738bfd3b3a4eeefd763dec15715d9d/thomaskuestner?format=bibtex'); return false;" href="https://puma.ub.uni-stuttgart.de/bibtex/2d4738bfd3b3a4eeefd763dec15715d9d/thomaskuestner?format=bibtex">[BibTeX]</a>&nbsp;<a onclick="toggleEndnote('thomaskuestner', 'd4738bfd3b3a4eeefd763dec15715d9d', 'https://puma.ub.uni-stuttgart.de/bibtex/2d4738bfd3b3a4eeefd763dec15715d9d/thomaskuestner?format=endnote'); return false;" href="https://puma.ub.uni-stuttgart.de/bibtex/2d4738bfd3b3a4eeefd763dec15715d9d/thomaskuestner?format=endnote">[Endnote]</a>

--------------------------------------------------------
Please read LICENSE file for licensing details.

Detailed information and installation instructions are available at: <br/>
https://sites.google.com/site/kspaceastronauts/compressed-sensing/cslab
