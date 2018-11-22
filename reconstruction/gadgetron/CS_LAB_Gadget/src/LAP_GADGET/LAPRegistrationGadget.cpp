/*
file name	: 	LAPRegistrationGadget.cpp
author		: 	Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)
version		: 	1.0
date		: 	15.01.2018
description	: 	implementation of the class LAPRegistrationGadget - only 4D data (x,y,z,t) provided
*/

#include "LAPRegistrationGadget.h"

#include "SomeFunctions.h"

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
	// create references for easier usage
	const hoNDArray<float> &data = *m2->getObjectPtr();
	std::vector<size_t> dimensions = *data.get_dimensions();

	// determine number of gate dimensions
	const size_t number_of_gate_dimensions = data.get_number_of_dimensions() > 3 ? data.get_number_of_dimensions()-3 : data.get_number_of_dimensions()-1;

	// get gate dimensions
	const std::vector<size_t> gate_dimensions(dimensions.end()-number_of_gate_dimensions, dimensions.end());

	// determine number of images
	const size_t number_of_images = std::accumulate(dimensions.end()-number_of_gate_dimensions, dimensions.end(), 1, std::multiplies<size_t>());

	// other dimensions are dimension of one image
	const std::vector<size_t> dimensions_of_image(dimensions.begin(), dimensions.end()-number_of_gate_dimensions);

	// check for image dimensions
	bool skip_registration = false;
	if (dimensions_of_image.size() < 2 || dimensions_of_image.size() > 3) {
		GWARN("Dataset does not contain 2D or 3D images - Skip registration...\n");
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

	// get fixed image from dataset
	std::vector<size_t> fixed_image_dimensions = dimensions_of_image;
	const hoNDArray<float> fixed_image(fixed_image_dimensions, const_cast<float*>(data.get_data_ptr()), false);

	// create registered (output) image
	hoNDArray<float> output_image(*data.get_dimensions());
	memcpy(output_image.get_data_ptr(), fixed_image.get_data_ptr(), fixed_image.get_number_of_bytes());

	// create array for deformation field with dimensions [X Y Z CombinedPhases VectorComponents(2 or 3)]
	std::vector<size_t> output_df_dimensions = dimensions_of_image;
	output_df_dimensions.push_back(number_of_images - 1);
	output_df_dimensions.push_back(dimensions_of_image.size());
	hoNDArray<float> output_deformation_field(&output_df_dimensions);

	CubeType cFixedImage = Cube<float>(dimensions_of_image.at(0), dimensions_of_image.at(1), dimensions_of_image.at(2));
	CubeType cMovingImage = Cube<float>(dimensions_of_image.at(0), dimensions_of_image.at(1), dimensions_of_image.at(2));

	memcpy(cFixedImage.memptr(), fixed_image.get_data_ptr(), fixed_image.get_number_of_bytes());

	//Construct the LocalAllpass Algorithm Object with Level min and max
	LAP3D mLAP3D(cFixedImage, cMovingImage, iLvlMin_, iLvlMax_);

	/* ------------------------------------------------------------------- */
	/* --------- loop over moving images and perform registration -------- */
	/* ------------------------------------------------------------------- */
	GINFO("Loop over moving images..\n");
	// loop over gates
// 	#pragma omp parallel for	// Parallelising here may work, but may also introduce errors. Check that before enabling!
	for (size_t calculated_images = 1; calculated_images < number_of_images; calculated_images++) {
		GINFO("%i of %i ...\n", calculated_images, number_of_images-1);

		// calculate offset for one single image
		const size_t moving_image_offset = std::accumulate(std::begin(dimensions_of_image), std::end(dimensions_of_image), 1, std::multiplies<size_t>())*calculated_images;  // accumulate() := dimensions_of_image[0]*dimensions_of_image[1]*...

		// crop moving image from 4D dataset
		std::vector<size_t> moving_image_dimensions = dimensions_of_image;
		hoNDArray<float> moving_image(moving_image_dimensions, const_cast<float*>(data.get_data_ptr()) + moving_image_offset, false);

		memcpy(cMovingImage.memptr(), moving_image.get_data_ptr(), moving_image.get_number_of_bytes());
		mLAP3D.setMovingImage(cMovingImage);

		// image registration
		field<CubeType> flow_estimation = mLAP3D.exec();

		// save deformation fields
		for (size_t i = 0; i < dimensions_of_image.size(); i++) {
			// calculate offset
			const size_t def_field_offset = moving_image.get_number_of_elements() * (	// all slots (4th dimension 3D images) have same size, so multiply it with position)
				i * output_deformation_field.get_size(3)	// select correct vector component (X|Y[|Z])
				+ (calculated_images - 1)	// select correct image
			);

			// copy image
			memcpy(output_deformation_field.get_data_ptr()+def_field_offset, flow_estimation(i).memptr(), moving_image.get_number_of_bytes());
		}

		// get output image
		// Shift first image according to estimated optical flow
		ShiftEngine3D shifter(cFixedImage, flow_estimation(0), flow_estimation(1), flow_estimation(2));
		CubeType cRegisteredImage = shifter.execCubicShift();

		// copy image to new registered 4D image
		memcpy(output_image.get_data_ptr()+moving_image_offset, cRegisteredImage.memptr(), moving_image.get_number_of_bytes());
	}

	// free memory
	m2->release();

	// new GadgetContainer
	GadgetContainerMessage<hoNDArray<float> > *cm2 = new GadgetContainerMessage<hoNDArray<float> >();

	// concatenate data with header
	m1->cont(cm2);

	// create output
	try {
		cm2->getObjectPtr()->create(*output_image.get_dimensions());
	} catch (std::runtime_error &err) {
		GEXCEPTION(err,"Unable to allocate new image array\n");
		m1->release();

		return GADGET_FAIL;
	}

	// copy data
	memcpy(cm2->getObjectPtr()->get_data_ptr(), output_image.begin(), output_image.get_number_of_bytes());

	// Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}

	// ########## create deformation field message ###########
	// copy header
	GadgetContainerMessage<ISMRMRD::ImageHeader> *m1_df = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
	fCopyImageHeader(m1_df, m1->getObjectPtr());

	// new GadgetContainer
	GadgetContainerMessage<hoNDArray<float> > *m2_df = new GadgetContainerMessage<hoNDArray<float> >();

	// concatenate data with header
	m1_df->cont(m2_df);

	// create output
	try {
		m2_df->getObjectPtr()->create(*output_deformation_field.get_dimensions());
	} catch (std::runtime_error &err) {
		GEXCEPTION(err,"Unable to allocate new deformation field array\n");
		m1_df->release();

		return GADGET_FAIL;
	}

	// copy data
	memcpy(m2_df->getObjectPtr()->get_data_ptr(), output_deformation_field.begin(), output_deformation_field.get_number_of_bytes());

	// Now pass on deformation field
	if (this->next()->putq(m1_df) < 0) {
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(LAPRegistrationGadget)
