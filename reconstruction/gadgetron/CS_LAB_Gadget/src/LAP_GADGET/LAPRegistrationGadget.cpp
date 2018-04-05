/*
file name	: 	LAPRegistrationGadget.cpp
author		: 	Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)
version		: 	1.0
date		: 	15.01.2018
description	: 	implementation of the class LAPRegistrationGadget - only 4D data (x,y,z,t) provided
*/

#include "LAPRegistrationGadget.h"

using namespace Gadgetron;

LAPRegistrationGadget::LAPRegistrationGadget()
{
}

LAPRegistrationGadget::~LAPRegistrationGadget()
{
}

int LAPRegistrationGadget::process_config(ACE_Message_Block *mb)
{
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	iLvlMin_ = LvlMin.value();
	iLvlMax_ = LvlMax.value();
#else
	iLvlMin_ = *(get_int_value("LvlMin").get());
	iLvlMax_ = *(get_int_value("LvlMax").get());
#endif

	return GADGET_OK;
}

int LAPRegistrationGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2)
{
	// create dimension variables
	// first image is fixed image (end-exhale position) all other images declared to be moving images
	std::vector<size_t> dimensions = *m2->getObjectPtr()->get_dimensions();

	// last dimension is number of images
	size_t number_of_images = dimensions.at(dimensions.size()-1);

	// other dimensions are dimension of one image
	std::vector<size_t> dimensions_of_image(dimensions.begin(), dimensions.end()-1);

	// check for image dimensions
	bool skip_registration = false;
	if (dimensions_of_image.size() != 3) {
		GWARN("Dataset does not contain a 3D image - Skip registration...\n");
		skip_registration = true;
	} else if (number_of_images <= 1) {
		GWARN("There must be at least to images to register them - Skip registration...\n");
		skip_registration = true;
	}

	// skip registration if it cannot be performed
	if (skip_registration) {
		if (this->next()->putq(m1) < 0) {
			return GADGET_FAIL;
		}

		return GADGET_OK;
	}

	// print dimensions
	std::stringstream ss;
	ss << "size - ";
	for (size_t i = 0; i < dimensions_of_image.size(); i++) {
		ss << dimensions_of_image.at(i) << " ";
	}
	ss << number_of_images;
	GDEBUG("%s\n", ss.str().c_str());

	// multiply all other elements in vector (will be number of pixels in one image)
	const size_t cuiNumberOfPixels = std::accumulate(std::begin(dimensions_of_image), std::end(dimensions_of_image), 1, std::multiplies<size_t>());

	float *pfDataset = m2->getObjectPtr()->get_data_ptr();

	// get fixed image from 4D dataset
	hoNDArray<float> fFixedImage(dimensions_of_image, pfDataset, false);

	// registered 4D image
	hoNDArray<float> fRegisteredImage(m2->getObjectPtr()->get_dimensions());
	memcpy(fRegisteredImage.get_data_ptr(), fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));

	CubeType cFixedImage = Cube<float>(dimensions_of_image.at(0), dimensions_of_image.at(1), dimensions_of_image.at(2));
	CubeType cMovingImage = Cube<float>(dimensions_of_image.at(0), dimensions_of_image.at(1), dimensions_of_image.at(2));

	memcpy(cFixedImage.memptr(), fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));

	//Construct the LocalAllpass Algorithm Object with Level min and max
	LAP3D mLAP3D(cFixedImage, cMovingImage, iLvlMin_, iLvlMax_);

	/* ------------------------------------------------------------------- */
	/* --------- loop over moving images and perform registration -------- */
	/* ------------------------------------------------------------------- */
	GINFO("Loop over moving images..\n");
	// loop over respiration
// 	#pragma omp parallel for	// Parallelising here may work, but may also introduce errors. Check that before enabling!
	for (size_t iState = 1; iState < number_of_images; iState++) {
		GINFO("%i of %i ...\n", iState, number_of_images-1);

		// crop moving image from 4D dataset
		size_t tOffset = std::accumulate(std::begin(dimensions_of_image), std::end(dimensions_of_image), 1, std::multiplies<size_t>())*iState;	// accumulate() := dimensions_of_image[0]*dimensions_of_image[1]*...
		hoNDArray<float> fMovingImage(dimensions_of_image, pfDataset + tOffset, false);

		memcpy(cMovingImage.memptr(), fMovingImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
		mLAP3D.setMovingImage(cMovingImage);

		// image registration
		field<CubeType> flow_estimation = mLAP3D.exec();

		// get output image
		// Shift first image according to estimated optical flow
		ShiftEngine3D shifter(cFixedImage, flow_estimation(0), flow_estimation(1), flow_estimation(2));
		CubeType cRegisteredImage = shifter.execCubicShift();

		// copy image to new registered 4D image
		memcpy(fRegisteredImage.get_data_ptr()+tOffset, cRegisteredImage.memptr(), cuiNumberOfPixels*sizeof(float));
	}

	// free memory
	m2->release();

	// new GadgetContainer
	GadgetContainerMessage<hoNDArray<float> > *cm2 = new GadgetContainerMessage<hoNDArray<float> >();

	// concatenate data with header
	m1->cont(cm2);

	// create output
	try {
		cm2->getObjectPtr()->create(*fRegisteredImage.get_dimensions());
	} catch (std::runtime_error &err) {
		GEXCEPTION(err,"Unable to allocate new image array\n");
		m1->release();

		return GADGET_FAIL;
	}

	// copy data
	memcpy(cm2->getObjectPtr()->get_data_ptr(), fRegisteredImage.begin(), sizeof(float)*fRegisteredImage.get_number_of_elements());

	// Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(LAPRegistrationGadget)
