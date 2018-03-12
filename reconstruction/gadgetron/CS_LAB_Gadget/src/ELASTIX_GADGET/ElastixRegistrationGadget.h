/*	
file name	: 	ElastixRegistrationGadget.h

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
				Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)

version		: 	1.2

date		: 	13.10.2015
				23.03.2017
				23.01.2018

description	: 	Elastix-based image registration
*/

#ifndef ELASTIXREGISTRATIONGADGET_H
#define ELASTIXREGISTRATIONGADGET_H

#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
	//#define GET_MACRO(_1,_2,_3,NAME,...) NAME
	//#define GDEBUG(...) GET_MACRO(__VA_ARGS__, GADGET_DEBUG1, GADGET_DEBUG2)(__VA_ARGS__)
	#define GADGET_DEBUG1(...) GDEBUG(__VA_ARGS__)
	#define GADGET_DEBUG2(x, ...) GDEBUG(x, ##__VA_ARGS__)
	#define GADGET_DEBUG_EXCEPTION(x,y) GEXCEPTION(x,y)
#endif

#pragma once
#include "CS_LAB_export.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include <ismrmrd.h>
#include "GadgetIsmrmrdReadWrite.h"

#include "elastixlib.h"				// elastix library
#include "itkParameterFileParser.h" // read parameter file from disk
#include "itkImage.h"				// itk image class
#include "itkImportImageFilter.h"	// itk image filter class

typedef elastix::ELASTIX::ParameterMapType RegistrationParametersType;
typedef itk::ParameterFileParser ParserType;
typedef itk::Image< float, 3 >  ImageType;
typedef itk::ImportImageFilter<float, 3> ImportFilterType;

namespace Gadgetron
{
	class EXPORTCSLAB ElastixRegistrationGadget : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<float>>
	{
		public:
		ElastixRegistrationGadget();
		~ElastixRegistrationGadget();
		
		GADGET_DECLARE(CS_Retro_NavigatorGadget);
		
		int process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2);
		int process_config(ACE_Message_Block* mb);		

		int fRegistration3D( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2);
		int fRegistration4D( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2);

		std::vector<size_t> vtDim_;
		bool bIs2D_;
		bool bIs3D_;
		bool bIs4D_;

		std::string sPathParam_;
		std::string sPathLog_;

		#if __GADGETRON_VERSION_HIGHER_3_6__ == 1
			GADGET_PROPERTY(PathParam, std::string, "PathParam", "");
			GADGET_PROPERTY(PathLog, std::string, "PathLog", "");
		#endif
	};
}
#endif //ELASTIXREGISTRATIONGADGET_H
