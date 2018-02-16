/*	
file name	: 	ElastixRegistrationGadget.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)
				Thomas Kuestner (thomas.kuestner@med.uni-tuebingen.de)

version		: 	1.2

date		: 	13.10.2015
				23.03.2017 - update for Gadgetron v3.8
				23.01.2018 - renaming

description	: 	implementation of the class ElastixRegistrationGadget - only 3D (x,y,t) and 4D data (x,y,z,t) provided

*/

#include "ElastixRegistrationGadget.h"

using namespace Gadgetron;

ElastixRegistrationGadget::ElastixRegistrationGadget():bIs2D_(false), bIs3D_(false), bIs4D_(false){};
ElastixRegistrationGadget::~ElastixRegistrationGadget(){};

int ElastixRegistrationGadget::process_config(ACE_Message_Block* mb){

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	sPathParam_ = PathParam.value();
	sPathLog_   = PathLog.value();
#else
	sPathParam_ = *(get_string_value("PathParam").get());
	sPathLog_   = *(get_string_value("PathLog").get());
#endif

	return GADGET_OK;
};

int ElastixRegistrationGadget::process( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2){
	
	/* ------------------------------------------------------------------- */
	/* --------------- check dimension of incoming dataset --------------- */
	/* ------------------------------------------------------------------- */
	std::vector<size_t> dimension = *m2->getObjectPtr()->get_dimensions();
	size_t num_dims = m2->getObjectPtr()->get_number_of_dimensions();
	size_t num_rep = 0, num_par = 0;
	
	// get dimensions flag
	if (num_dims == 2){
		bIs2D_ = true;
	}
	else if (num_dims == 3){
		bIs3D_ = true;
	}
	else if (num_dims == 4){
		bIs4D_ = true;
	}	

	/* ------------------------------------------------------------------- */
	/* ----------------------- call registration ------------------------- */
	/* ------------------------------------------------------------------- */
	if (bIs2D_){		
		GWARN("2D dataset detected..unable to perform registration - step skipped!\n");
	}
	else if(bIs3D_){
		GINFO("3D dataset detected..perform registration in 3rd dimension!\n");
		fRegistration3D(m1, m2);
	}
	else if(bIs4D_){
		GINFO("4D dataset detected..perform registration in 4th dimension!\n");
		fRegistration4D(m1, m2);
	}
	else{
		GWARN("unknown dataset detected..unable to perform registration - step skipped!\n");
	}

	//Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}
	
	return GADGET_OK;
};

int ElastixRegistrationGadget::fRegistration3D( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2){
	/* ------------------------------------------------------------------- */
	/* --------------- create registration parameter data ---------------- */
	/* ------------------------------------------------------------------- */
	GINFO("Load elastix parameter file..\n");
	// Create parser for transform parameters text file.
	ParserType::Pointer file_parser = ParserType::New();
	// Try parsing transform parameters text file.
	GINFO("search for parameter file - %s..\n", sPathParam_.c_str());
	file_parser->SetParameterFileName( sPathParam_ );
	try {
		file_parser->ReadParameterFile();
	}
	catch( itk::ExceptionObject & e ) {
		std::cout << e.what() << std::endl;
	};
	RegistrationParametersType parameters = file_parser->GetParameterMap();
	GDEBUG("parameter file - %s - loaded..\n", sPathParam_.c_str());

	/* ------------------------------------------------------------------- */
	/* -------------------- init elastix registration -------------------- */
	/* ------------------------------------------------------------------- */
	// first image is fixed image all other images declared to be moving images
	vtDim_ = *m2->getObjectPtr()->get_dimensions();

	int iNoImages = vtDim_[2];
	const unsigned int cuiNumberOfPixels = vtDim_[0]*vtDim_[1];
	float *pfDataset = m2->getObjectPtr()->get_data_ptr();
	
	// get fixed image from 3D dataset
	hoNDArray<float> fFixedImage(vtDim_[0], vtDim_[1], pfDataset, false);
		
	// registered 3D image
	hoNDArray<float> fRegisteredImage(m2->getObjectPtr()->get_dimensions());
	memcpy(fRegisteredImage.get_data_ptr(), fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));

	// set itk image parameter
	ImportFilterType::Pointer itkImportFilter = ImportFilterType::New();
	ImportFilterType::SizeType itkSize;
	itkSize[0] = vtDim_[0];
	itkSize[1] = vtDim_[1];
	itkSize[2] = 1;
	
	ImportFilterType::IndexType itkStart;
	itkStart.Fill(0);

	ImportFilterType::RegionType itkRegion;
	itkRegion.SetIndex(itkStart);
	itkRegion.SetSize(itkSize);
		
	itkImportFilter->SetRegion(itkRegion);
	
	double itkOrigin[3];
	itkOrigin[0] = 0.0;
	itkOrigin[1] = 0.0;
	itkOrigin[2] = 0.0;
	itkImportFilter->SetOrigin(itkOrigin);

	double itkSpacing[3];
	itkSpacing[0] = m1->getObjectPtr()->field_of_view[0]/m1->getObjectPtr()->matrix_size[0];
	itkSpacing[1] = m1->getObjectPtr()->field_of_view[1]/m1->getObjectPtr()->matrix_size[1];
	itkSpacing[2] = m1->getObjectPtr()->field_of_view[2]/m1->getObjectPtr()->matrix_size[2];
	itkImportFilter->SetSpacing(itkSpacing);
	
	float *fLocalBuffer = new float[cuiNumberOfPixels];
	memcpy(fLocalBuffer, fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
	itkImportFilter->SetImportPointer(fLocalBuffer, cuiNumberOfPixels, true);	
	ImageType::Pointer itkFixedImage = ImageType::New();
	itkFixedImage = itkImportFilter->GetOutput();
	itkImportFilter->Update();
	
	/* ------------------------------------------------------------------- */
	/* --------- loop over moving images and perform registration -------- */
	/* ------------------------------------------------------------------- */
	GINFO("Registration of images..\n");
	
	// loop over respiration 
	for (int iState = 1; iState < iNoImages; iState++){
		GINFO("%i of %i ...\n", iState, iNoImages-1);
		// crop moving image from 4D dataset
		size_t tOffset = vtDim_[0]*vtDim_[1]*iState;
		hoNDArray<float> fMovingImage(vtDim_[0], vtDim_[1], pfDataset + tOffset, false);
		
		float *fLocalBuffer = new float[cuiNumberOfPixels];
		memcpy(fLocalBuffer, fMovingImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
		itkImportFilter->SetImportPointer(fLocalBuffer, cuiNumberOfPixels, true);	
		ImageType::Pointer itkMovingImage = ImageType::New();
		itkMovingImage = itkImportFilter->GetOutput();
		itkImportFilter->Update();
		
		// image registration
		elastix::ELASTIX* elastix_obj = new elastix::ELASTIX();
		int error = 0;
		try{
			// perform registration with (FixedImage, MovingImage, ParameterFile, OutputPath, ElastixLog, ConsoleOutput, FixedImageMask, MovingImageMask)
			error = elastix_obj->RegisterImages(static_cast<itk::DataObject::Pointer>(itkFixedImage.GetPointer()), static_cast<itk::DataObject::Pointer>(itkMovingImage.GetPointer()), parameters, sPathLog_, true, false, 0, 0);
		}
		catch (itk::ExitEvent &err){
			// error handling - write message and fill array with zeros
			GERROR("Error event catched directly from elastix\n");
		}
		
		// get output image
		float *fLocalBufferOutput = new float[cuiNumberOfPixels];
		ImageType * itkOutputImage;
		if( error == 0 ) {
			if( elastix_obj->GetResultImage().IsNotNull() ){
				itkOutputImage = static_cast<ImageType*>(elastix_obj->GetResultImage().GetPointer());				
			}
			else{
				GERROR("GetResultImage() is NULL \n", error);
			}
		}
		else {
			// error handling - write message and fill array with zeros
			GERROR("array is zero\n", error);
			itkOutputImage->FillBuffer(0.0);
		}
		
		// copy image to new registered 3D image
		memcpy(fRegisteredImage.get_data_ptr()+tOffset, itkOutputImage->GetBufferPointer(), cuiNumberOfPixels*sizeof(float));

		// clean up
		delete elastix_obj;
	}
	
	// new GadgetContainer
	GadgetContainerMessage< hoNDArray< float > >* cm2 = new GadgetContainerMessage<hoNDArray< float > >();
    
	// concatenate data with header
	m1->cont(cm2);
	
	// create output
	try{cm2->getObjectPtr()->create(&vtDim_);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to allocate new image array\n");
		m1->release();
		return -1;
	}
	
	// copy data
	memcpy(cm2->getObjectPtr()->get_data_ptr(), fRegisteredImage.begin(), sizeof(float)*fRegisteredImage.get_number_of_elements());

	return GADGET_OK;
};

int ElastixRegistrationGadget::fRegistration4D( GadgetContainerMessage< ISMRMRD::ImageHeader>* m1, GadgetContainerMessage< hoNDArray< float > >* m2){
	/* ------------------------------------------------------------------- */
	/* --------------- create registration parameter data ---------------- */
	/* ------------------------------------------------------------------- */
	GINFO("Load elastix parameter file..\n");
	// Create parser for transform parameters text file.
	ParserType::Pointer file_parser = ParserType::New();
	// Try parsing transform parameters text file.
	GINFO("search for parameter file - %s..\n", sPathParam_.c_str());
	file_parser->SetParameterFileName( sPathParam_ );
	try {
		file_parser->ReadParameterFile();
	}
	catch( itk::ExceptionObject & e ) {
		std::cout << e.what() << std::endl;
	};
	RegistrationParametersType parameters = file_parser->GetParameterMap();
	GINFO("parameter file - %s - loaded..\n", sPathParam_.c_str());

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

	// set itk image parameter
	ImportFilterType::Pointer itkImportFilter = ImportFilterType::New();
	ImportFilterType::SizeType itkSize;
	itkSize[0] = vtDim_[0];
	itkSize[1] = vtDim_[1];
	itkSize[2] = vtDim_[2];
	
	ImportFilterType::IndexType itkStart;
	itkStart.Fill(0);

	ImportFilterType::RegionType itkRegion;
	itkRegion.SetIndex(itkStart);
	itkRegion.SetSize(itkSize);
		
	itkImportFilter->SetRegion(itkRegion);
	
	double itkOrigin[3];
	itkOrigin[0] = 0.0;
	itkOrigin[1] = 0.0;
	itkOrigin[2] = 0.0;
	itkImportFilter->SetOrigin(itkOrigin);

	double itkSpacing[3];
	itkSpacing[0] = m1->getObjectPtr()->field_of_view[0]/m1->getObjectPtr()->matrix_size[0];
	itkSpacing[1] = m1->getObjectPtr()->field_of_view[1]/m1->getObjectPtr()->matrix_size[1];
	itkSpacing[2] = m1->getObjectPtr()->field_of_view[2]/m1->getObjectPtr()->matrix_size[2];
	itkImportFilter->SetSpacing(itkSpacing);
	
	float *fLocalBuffer = new float[cuiNumberOfPixels];
	memcpy(fLocalBuffer, fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
	itkImportFilter->SetImportPointer(fLocalBuffer, cuiNumberOfPixels, true);	
	ImageType::Pointer itkFixedImage = ImageType::New();
	itkFixedImage = itkImportFilter->GetOutput();
	itkImportFilter->Update();
	
	/* ------------------------------------------------------------------- */
	/* --------- loop over moving images and perform registration -------- */
	/* ------------------------------------------------------------------- */
	GINFO("Loop over moving images..\n");
	// loop over respiration 
	for (int iState = 1; iState < iNoImages; iState++){
		GINFO("%i of %i ...\n", iState, iNoImages-1);
		// crop moving image from 4D dataset
		size_t tOffset = vtDim_[0]*vtDim_[1]*vtDim_[2]*iState;
		hoNDArray<float> fMovingImage(vtDim_[0], vtDim_[1], vtDim_[2], pfDataset + tOffset, false);
		
		float *fLocalBuffer = new float[cuiNumberOfPixels];
		memcpy(fLocalBuffer, fMovingImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
		itkImportFilter->SetImportPointer(fLocalBuffer, cuiNumberOfPixels, true);	
		ImageType::Pointer itkMovingImage = ImageType::New();
		itkMovingImage = itkImportFilter->GetOutput();
		itkImportFilter->Update();
		
		// image registration
		elastix::ELASTIX* elastix_obj = new elastix::ELASTIX();
		int error = 0;
		try{
			// perform registration with (FixedImage, MovingImage, ParameterFile, OutputPath, ElastixLog, ConsoleOutput, FixedImageMask, MovingImageMask)
			error = elastix_obj->RegisterImages(static_cast<itk::DataObject::Pointer>(itkFixedImage.GetPointer()), static_cast<itk::DataObject::Pointer>(itkMovingImage.GetPointer()), parameters, "D:/Programme/elastix", true, false, 0, 0);
		}
		catch (itk::ExitEvent &err){
			// error handling - write message and fill array with zeros
			GERROR("Error event catched directly from elastix\n");
		}
		
		// get output image
		float *fLocalBufferOutput = new float[cuiNumberOfPixels];
		ImageType * itkOutputImage;
		if( error == 0 ) {
			if( elastix_obj->GetResultImage().IsNotNull() ){
				itkOutputImage = static_cast<ImageType*>(elastix_obj->GetResultImage().GetPointer());				
			}
			else{
				GERROR("GetResultImage() is NULL \n", error);
			}
		}
		else {
			// error handling - write message and fill array with zeros
			GERROR("array is zero\n", error);
			itkOutputImage->FillBuffer(0.0);
		}
		
		// copy image to new registered 4D image
		memcpy(fRegisteredImage.get_data_ptr()+tOffset, itkOutputImage->GetBufferPointer(), cuiNumberOfPixels*sizeof(float));

		// clean up
		delete elastix_obj;
	}
	
	// new GadgetContainer
	GadgetContainerMessage< hoNDArray< float > >* cm2 = new GadgetContainerMessage<hoNDArray< float > >();
    
	// concatenate data with header
	m1->cont(cm2);
	
	// create output
	try{cm2->getObjectPtr()->create(&vtDim_);}
	catch (std::runtime_error &err){
		GEXCEPTION(err,"Unable to allocate new image array\n");
		m1->release();
		return -1;
	}
	
	// copy data
	memcpy(cm2->getObjectPtr()->get_data_ptr(), fRegisteredImage.begin(), sizeof(float)*fRegisteredImage.get_number_of_elements());

	return GADGET_OK;
};

GADGET_FACTORY_DECLARE(ElastixRegistrationGadget)
