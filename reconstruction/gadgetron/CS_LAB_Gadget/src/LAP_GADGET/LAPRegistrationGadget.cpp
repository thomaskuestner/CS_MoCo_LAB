/*	
file name	: 	LAPRegistrationGadget.cpp

author		: 	Thomas Kuestner	(thomas.kuestner@med.uni-tuebingen.de)

version		: 	1.0

date		: 	15.01.2018

description	: 	implementation of the class LAPRegistrationGadget - only 4D data (x,y,z,t) provided

*/

#include "LAPRegistrationGadget.h"

using namespace Gadgetron;

LAPRegistrationGadget::LAPRegistrationGadget():bIs2D_(false), bIs3D_(false), bIs4D_(false), iLvlMin_(0), iLvlMax_(4)
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
	/* ------------------------------------------------------------------- */
	/* --------------- check dimension of incoming dataset --------------- */
	/* ------------------------------------------------------------------- */
	std::vector<size_t> dimension = *m2->getObjectPtr()->get_dimensions();
	size_t num_dims = m2->getObjectPtr()->get_number_of_dimensions();

	// get dimensions flag
	if (num_dims == 2) {
		bIs2D_ = true;
	} else if (num_dims == 3) {
		bIs3D_ = true;
	} else if (num_dims == 4) {
		bIs4D_ = true;
	}

	/* ------------------------------------------------------------------- */
	/* ----------------------- call registration ------------------------- */
	/* ------------------------------------------------------------------- */
	if (bIs2D_) {
		GWARN("2D dataset detected..unable to perform registration - step skipped!\n");
	} else if (bIs3D_) {
		GWARN("3D dataset detected..unable to perform registration - step skipped!\n");
	} else if (bIs4D_) {
		GINFO("4D dataset detected..perform registration in 4th dimension!\n");
		fRegistration4D(m1, m2);
	} else {
		GWARN("unknown dataset detected..unable to perform registration - step skipped!\n");
	}

	//Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}
	
	return GADGET_OK;
};

int LAPRegistrationGadget::fRegistration4D(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2)
{
	// first image is fixed image (end-exhale position) all other images declared to be moving images
	vtDim_ = *m2->getObjectPtr()->get_dimensions();

	GDEBUG("size - %i %i %i %i\n", vtDim_[0], vtDim_[1], vtDim_[2], vtDim_[3]);

	int iNoImages = vtDim_[3];
	const unsigned int cuiNumberOfPixels = vtDim_[0]*vtDim_[1]*vtDim_[2];
	float *pfDataset = m2->getObjectPtr()->get_data_ptr();

	// get fixed image from 4D dataset
	hoNDArray<float> fFixedImage(vtDim_[0], vtDim_[1], vtDim_[2], pfDataset, false);

	// registered 4D image
	hoNDArray<float> fRegisteredImage(m2->getObjectPtr()->get_dimensions());
	memcpy(fRegisteredImage.get_data_ptr(), fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));

	CubeType cFixedImage = Cube<float>(vtDim_[0], vtDim_[1], vtDim_[2]);
	CubeType cMovingImage = Cube<float>(vtDim_[0], vtDim_[1], vtDim_[2]);

	memcpy(cFixedImage.memptr(), fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
	
	//Construct the LocalAllpass Algorithm Object with Level min and max
	LAP3D mLAP3D(cFixedImage, cMovingImage, iLvlMin_, iLvlMax_);

	/* ------------------------------------------------------------------- */
	/* --------- loop over moving images and perform registration -------- */
	/* ------------------------------------------------------------------- */
	GINFO("Loop over moving images..\n");
	// loop over respiration 
	for (int iState = 1; iState < iNoImages; iState++) {
		GINFO("%i of %i ...\n", iState, iNoImages-1);

		// crop moving image from 4D dataset
		size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iState;
		hoNDArray<float> fMovingImage(vtDim_[0], vtDim_[1], vtDim_[2], pfDataset + tOffset, false);

		memcpy(cMovingImage.memptr(), fMovingImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
		mLAP3D.setMovingImage(cMovingImage);

		// image registration
		field<CubeType> flow_estimation = mLAP3D.exec();

		// get output image
		//Shift first image according to estimated optical flow
		ShiftEngine3D shifter(cFixedImage, flow_estimation(0), flow_estimation(1), flow_estimation(2));
		CubeType cRegisteredImage = shifter.execCubicShift();

		// copy image to new registered 4D image
		memcpy(fRegisteredImage.get_data_ptr()+tOffset, cRegisteredImage.memptr(), cuiNumberOfPixels*sizeof(float));
	}

	// new GadgetContainer
	GadgetContainerMessage<hoNDArray<float> > *cm2 = new GadgetContainerMessage<hoNDArray<float> >();

	// concatenate data with header
	m1->cont(cm2);

	// create output
	try {
		cm2->getObjectPtr()->create(&vtDim_);
	} catch (std::runtime_error &err) {
		GEXCEPTION(err,"Unable to allocate new image array\n");
		m1->release();
		return -1;
	}

	// copy data
	memcpy(cm2->getObjectPtr()->get_data_ptr(), fRegisteredImage.begin(), sizeof(float)*fRegisteredImage.get_number_of_elements());

	return GADGET_OK;
};

GADGET_FACTORY_DECLARE(LAPRegistrationGadget)
