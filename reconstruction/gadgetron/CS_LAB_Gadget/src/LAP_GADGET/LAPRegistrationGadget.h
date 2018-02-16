/*	
file name	: 	LAPRegistrationGadget.h

author		: 	Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)

version		: 	1.0

date		: 	15.01.2018

description	: 	LAP-based image registration
*/

#ifndef LAPREGISTRATIONGADGET_H
#define LAPREGISTRATIONGADGET_H

#include "gadgetron_messages.h"

#pragma once
#include "CS_LAB_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include <ismrmrd.h>
#include "GadgetIsmrmrdReadWrite.h"

#include "lap3d.h"
#include "shiftengine3d.h"

#include "itkParameterFileParser.h" // read parameter file from disk
#include "itkImage.h"				// itk image class
//#include "itkImportImageFilter.h"	// itk image filter class

typedef itk::ParameterFileParser ParserType;
typedef itk::Image< float, 3 >  ImageType;
//typedef itk::ImportImageFilter<float, 3> ImportFilterType;

namespace Gadgetron
{
	class EXPORTCSLAB LAPRegistrationGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<float>>
	{
		public:
		LAPRegistrationGadget();
		~LAPRegistrationGadget();
		
		GADGET_DECLARE(LAPRegistrationGadget);
		
		int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2);
		int process_config(ACE_Message_Block* mb);		

		int fRegistration4D( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2);

		std::vector<size_t> vtDim_;
		bool bIs2D_;
		bool bIs3D_;
		bool bIs4D_;

		int iLvlMin_;
		int iLvlMax_;

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		GADGET_PROPERTY(LvlMin, int, "LvlMin", 0);
		GADGET_PROPERTY(LvlMax, int, "LvlMax", 4);
#endif
	};
}
#endif //LAPREGISTRATIONGADGET_H
