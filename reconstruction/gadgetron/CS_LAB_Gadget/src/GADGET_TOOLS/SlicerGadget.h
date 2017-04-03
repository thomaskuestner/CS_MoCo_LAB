/*	
file name	: 	SlicerGadget.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	cuts the 3D/4D data into 2D images before sending them to the ICE pipeline. For a detailed description it is suggested to read p. 32-34 of the thesis.

input		:	mb					: 	ACE_Message_Block, which includes the flexible header data
				m1					: 	appended image header information
				m2					: 	2D, 2Dt, 3D, 4D input image data
				
output		:	cm1_sec_buffer_		: 	new header information
				sec_buffer_			: 	new 2D image data
				
functions	:	SlicerGadget()		: 	class constructor
				~SlicerGadget()		: 	class destructor
				process_config(...)	: 	reads the flexible header data
				process(...)		:	data processing - cuts the data into 2D data sets
				GADGET_DECLARE(...)	:	Gadget declaration (for the Gadgetron pipeline)
				
variables	:	image_counter_		:	running image index
				image_series_		:	running series index
				bIs2D_				:	flag, which indicates 2D data set
				bIs3D_				:	flag, which indicates 3D data set
				bIs4D_				:	flag, which indicates 4D data set
				
references	:	-
*/
#ifndef SLICERGADGET_H
#define SLICERGADGET_H

#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "CS_LAB_export.h"
#include <ismrmrd.h>
#include <complex>
#include "GadgetIsmrmrdReadWrite.h"
#include "GlobalVar.h"
#include "SomeFunctions.h"

namespace Gadgetron{
  
	class EXPORTCSLAB SlicerGadget : public Gadget2< ISMRMRD::ImageHeader, hoNDArray<float> >
    {
    public:      
      SlicerGadget();
      ~SlicerGadget();
      int process_config(ACE_Message_Block* mb);
	  int process(GadgetContainerMessage<ISMRMRD::ImageHeader>*m1, GadgetContainerMessage<hoNDArray<float>>*m2);
	  GADGET_DECLARE(SlicerGadget);

    protected:
      long long image_counter_;
      long long image_series_;
	  bool bIs2D_;
	  bool bIs3D_;
	  bool bIs4D_;
    };
}
#endif //SLICERGADGET_H
