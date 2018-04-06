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

#include <boost/algorithm/string/predicate.hpp>
#include "SomeFunctions.h"

using namespace Gadgetron;

ElastixRegistrationGadget::ElastixRegistrationGadget()
{
}

ElastixRegistrationGadget::~ElastixRegistrationGadget()
{
}

int ElastixRegistrationGadget::process_config(ACE_Message_Block *mb)
{
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	log_output_	= LogOutput.value();
	sPathParam_	= PathParam.value();
	sPathLog_	= PathLog.value();
#else
	log_output_	= *(get_bool_value("LogOutput").get());
	sPathParam_	= *(get_string_value("PathParam").get());
	sPathLog_	= *(get_string_value("PathLog").get());
#endif

	// ensure that path variable ends with dir separator
#ifdef __WIN32__
	#define DIR_SEPARATOR "\\"
#else
	#define DIR_SEPARATOR "/"
#endif

	if (boost::algorithm::ends_with(sPathLog_, DIR_SEPARATOR)) {
		sPathLog_ += DIR_SEPARATOR;
	}

	return GADGET_OK;
}

int ElastixRegistrationGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<float> > *m2)
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

	/* ------------------------------------------------------------------- */
	/* --------------- create registration parameter data ---------------- */
	/* ------------------------------------------------------------------- */
	GINFO("Load elastix parameter file..\n");

	// Create parser for transform parameters text file.
	ParserType::Pointer file_parser = ParserType::New();

	GINFO("search for parameter file - %s..\n", sPathParam_.c_str());

	// Try parsing transform parameters text file.
	file_parser->SetParameterFileName(sPathParam_);

	try {
		file_parser->ReadParameterFile();
	} catch (itk::ExceptionObject &e) {
		std::cout << e.what() << std::endl;
	}

	RegistrationParametersType parameters = file_parser->GetParameterMap();

	GINFO("parameter file - %s - loaded..\n", sPathParam_.c_str());

	/* ------------------------------------------------------------------- */
	/* -------------------- init elastix registration -------------------- */
	/* ------------------------------------------------------------------- */

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

	// get fixed image from dataset
	hoNDArray<float> fFixedImage(dimensions_of_image, pfDataset, false);

	// and create new registered image
	std::vector<size_t> new_reg_image_dim = *m2->getObjectPtr()->get_dimensions();
	new_reg_image_dim.at(new_reg_image_dim.size()-1) += (number_of_images - 1) * dimensions_of_image.size();	// append images for deformation fields (images - fixedImage) * dimensions per image (2D/3D)
	hoNDArray<float> fRegisteredImage(&new_reg_image_dim);
	memcpy(fRegisteredImage.get_data_ptr(), fFixedImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));

	// set itk image parameter
	ImportFilterType::Pointer itkImportFilter = ImportFilterType::New();
	ImportFilterType::SizeType itkSize;
	itkSize[0] = dimensions_of_image.at(0);
	itkSize[1] = dimensions_of_image.at(1);
	itkSize[2] = dimensions_of_image.size() >= 3 ? dimensions_of_image.at(2) : 1;	// ITK handels 3D image, so if we have a 2D moving image, 3rd dimension=1

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
// 	#pragma omp parallel for	// Note: Elastix itselfs parallels quite good. Some serious error can occur if you parallelise here!
	for (size_t iState = 1; iState < number_of_images; iState++) {
		GINFO("%i of %i ...\n", iState, number_of_images-1);

		// crop moving image from 4D dataset
		size_t tOffset = std::accumulate(std::begin(dimensions_of_image), std::end(dimensions_of_image), 1, std::multiplies<size_t>())*iState;	// accumulate() := dimensions_of_image[0]*dimensions_of_image[1]*...
		hoNDArray<float> fMovingImage(dimensions_of_image, pfDataset + tOffset, false);

		float *fLocalBuffer = new float[cuiNumberOfPixels];
		memcpy(fLocalBuffer, fMovingImage.get_data_ptr(), cuiNumberOfPixels*sizeof(float));
		itkImportFilter->SetImportPointer(fLocalBuffer, cuiNumberOfPixels, true);
		ImageType::Pointer itkMovingImage = ImageType::New();
		itkMovingImage = itkImportFilter->GetOutput();
		itkImportFilter->Update();

		// image registration
		elastix::ELASTIX *elastix_obj = new elastix::ELASTIX();
		int error = 0;

		try {
			// perform registration with (FixedImage, MovingImage, ParameterFile, OutputPath, ElastixLog, ConsoleOutput, FixedImageMask, MovingImageMask)
			error = elastix_obj->RegisterImages(static_cast<itk::DataObject::Pointer>(itkFixedImage.GetPointer()), static_cast<itk::DataObject::Pointer>(itkMovingImage.GetPointer()), parameters, sPathLog_, log_output_, false, 0, 0);
		} catch (itk::ExitEvent &err) {
			// error handling - write message and fill array with zeros
			GERROR("Error event catched directly from elastix\n");
		}

		// get output image
		ImageType *itkOutputImage = NULL;

		if (error == 0) {
			if (elastix_obj->GetResultImage().IsNotNull()) {
				itkOutputImage = static_cast<ImageType*>(elastix_obj->GetResultImage().GetPointer());
			} else {
				GERROR("GetResultImage() is NULL \n", error);
			}
		} else {
			// error handling - write message and fill array with zeros
			GERROR("array is zero\n", error);
		}

		// set output image to zeros if non-existent
		void *memcpy_pointer = NULL;
		ImageType::Pointer zero_image;
		if (itkOutputImage != NULL) {
			memcpy_pointer = itkOutputImage->GetBufferPointer();
		} else {
			ImageType::IndexType start;
			for (size_t i = 0; i < 3; i++) {
				start[i] = 0;
			}

			ImageType::SizeType size;
			for (size_t i = 0; i < dimensions_of_image.size(); i++) {
				size[i] = dimensions_of_image[i];
			}
			size[dimensions_of_image.size()] = number_of_images;

			ImageType::RegionType region;
			region.SetSize(size);
			region.SetIndex(start);

			zero_image = ImageType::New();
			zero_image->SetRegions(region);
			zero_image->Allocate();
			zero_image->FillBuffer(0.0);

			memcpy_pointer = zero_image->GetBufferPointer();
		}

		// copy image to new registered image
		memcpy(fRegisteredImage.get_data_ptr()+tOffset, memcpy_pointer, cuiNumberOfPixels*sizeof(float));

		// clean up
		delete elastix_obj;

		// calculate deformation field
		std::string transformix_command = std::string("transformix -def all -out ")+sPathLog_+std::string(" -tp ")+sPathLog_+std::string("TransformParameters.0.txt");
		int term_status = system(transformix_command.c_str());

		if (WIFEXITED(term_status)) {
			// transformix completed successfully. Handle case here.
			// create new dimension vector for array (dimensions e.g. 3 256 256 72) - should also apply on 2D image
			std::vector<size_t> deformation_field_dims;
			deformation_field_dims.push_back(dimensions_of_image.size());
			for (size_t i = 0; i < dimensions_of_image.size(); i++) {
				deformation_field_dims.push_back(dimensions_of_image.at(i));
			}

			// create new array for deformation field image
			hoNDArray<float> deformation_field(&deformation_field_dims);

			std::string deformation_field_file_name = sPathLog_ + std::string("deformationField.mhd");

			// read out image
			switch (dimensions_of_image.size()) {
			case 2:
				read_itk_to_hondarray<2>(deformation_field, deformation_field_file_name.c_str(), cuiNumberOfPixels);
				break;
			case 3:
				read_itk_to_hondarray<3>(deformation_field, deformation_field_file_name.c_str(), cuiNumberOfPixels);
				break;
			default:
				GERROR("Image dimension (%d) is neither 2D nor 3D! Abort.\n", dimensions_of_image.size());
				return GADGET_FAIL;
			}

			// permute data (in ITK dim is: [value_vector size_x size_y size_z])
			// we want to output: (layer with value_vector X component, layer with Y,...)
			// so: [value X Y Z] -> [X Y Z value]
			std::vector<size_t> new_dim_order;
			new_dim_order.push_back(1);
			new_dim_order.push_back(2);
			new_dim_order.push_back(3);
			new_dim_order.push_back(0);
			deformation_field = *permute(&deformation_field, &new_dim_order, false);

			// now calculate offset for the great return image
			size_t def_field_offset = cuiNumberOfPixels * (	// all slots (4th dimension 3D images) have same size, so multiply it with position)
				number_of_images							// jump over all valid images
				+ (iState - 1) * dimensions_of_image.size()	// jump to correct deformation field image position (*dim.size() because each DF Image consists of N parts (x,y[,z]))
			);

			// and copy permuted data into great return image
			memcpy(fRegisteredImage.get_data_ptr()+def_field_offset, deformation_field.get_data_ptr(), deformation_field.get_number_of_elements()*sizeof(float));
		} else {
			GWARN("No deformation field for state %d available (Errors in transformix).\n", iState);
		}
	}

	// clean up the elastix/transformix generated files
	std::vector<std::string> files_to_remove;
	files_to_remove.push_back(std::string("deformationField.mhd"));
	files_to_remove.push_back(std::string("deformationField.raw"));
	files_to_remove.push_back(std::string("transformix.log"));
	files_to_remove.push_back(std::string("TransformParameters.0.txt"));
	files_to_remove.push_back(std::string("IterationInfo.0.R0.txt"));
	files_to_remove.push_back(std::string("IterationInfo.0.R1.txt"));
	files_to_remove.push_back(std::string("IterationInfo.0.R2.txt"));
	files_to_remove.push_back(std::string("IterationInfo.0.R3.txt"));

	while (files_to_remove.size() > 0) {
		std::string file_to_remove = sPathLog_ + files_to_remove.back();

		if (remove(file_to_remove.c_str()) != 0) {
			GWARN("Could not remove %s. Please delete it manually!\n", file_to_remove.c_str());
		}

		// delete last element in list
		files_to_remove.pop_back();
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

	// correct channels value for MRIImageWriter (last dimension of output array)
	m1->getObjectPtr()->channels = fRegisteredImage.get_size(fRegisteredImage.get_number_of_dimensions()-1);

	// Now pass on image
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(ElastixRegistrationGadget)
