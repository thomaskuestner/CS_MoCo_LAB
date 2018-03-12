#include "ImageSaveHDFGadgetCPLX.h"

using namespace Gadgetron;

// process config - get filename, ..
int ImageSaveHDFGadgetCPLX::process_config(ACE_Message_Block *mb)
{
	// ace strings for debug output
	ACE_TCHAR ace_file_path[4096];
	ACE_TCHAR ace_ismrmrd_file_name[4096];

	file_prefix_ = *(get_string_value("file_prefix").get());
	if (file_prefix_.empty()) {
		file_prefix_ = "ISMRMRD_SAVE";
	}

	// get file path - if empty -> environment variable of gadgetron
	file_path_ = *(get_string_value("file_path").get());
	if (file_path_.empty()) {
		file_path_ = (std::string)getenv("GADGETRON_HOME");
		ACE_OS_String::strncpy(ace_file_path, file_path_.c_str(), 4096);
		ACE_DEBUG((LM_INFO, ACE_TEXT("\n -- no directory in xml config file - assume default path - file path changed to: %s\n"), ace_file_path));
	}

#ifdef __WIN32__
	if (!dirExists(file_path_)) {
		// create directory
		if (!CreateDirectory((LPCSTR)file_path_.c_str(), NULL)) {
			// set default value
			file_path_ = (std::string)getenv("GADGETRON_HOME");
			ACE_OS_String::strncpy(ace_file_path, file_path_.c_str(), 4096);
			ACE_DEBUG((LM_INFO, ACE_TEXT("\n -- failed to create directory - file path changed to: %s\n"), ace_file_path));
		}
	}
#else
	if (!stat(file_path_)) {
		// create directory
		if (mkdir(file_path_.c_str())) {
			// set default value
			file_path_ = (std::string)getenv("GADGETRON_HOME");
			ACE_OS_String::strncpy(ace_file_path, file_path_.c_str(), 4096);
			ACE_DEBUG((LM_INFO, ACE_TEXT("\n -- failed to create directory - file path changed to: %s\n"), ace_file_path));
		}
	}
#endif

	//Generate filename
	ismrmrd_file_name_ = file_path_ + file_prefix_ + std::string("_") + get_date_time_string() + std::string(".h5");
	ISMRMRD::HDF5Exclusive lock; //This will ensure threadsafe access to HDF5
	ismrmrd_dataset_ = boost::shared_ptr<ISMRMRD::IsmrmrdDataset>(new ISMRMRD::IsmrmrdDataset(ismrmrd_file_name_.c_str(), "dataset"));

	std::string xml_config(mb->rd_ptr());

	if (ismrmrd_dataset_->writeHeader(xml_config) < 0) {
		GERROR("Failed to write XML header to HDF file\n");

		return GADGET_FAIL;
	}

	ACE_OS_String::strncpy(ace_file_path, file_path_.c_str(), 4096);
	ACE_OS_String::strncpy(ace_ismrmrd_file_name, ismrmrd_file_name_.c_str(), 4096);
	ACE_DEBUG((LM_INFO, ACE_TEXT("\n -- file path         : %s\n"), ace_file_path));
	ACE_DEBUG((LM_INFO,   ACE_TEXT(" -- hdf5 file name    : %s\n"), ace_ismrmrd_file_name));

	return GADGET_OK;
}

// save the incoming acquisition data to an hdf5 file
int ImageSaveHDFGadgetCPLX::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	std::vector<size_t> dim = *m2->getObjectPtr()->get_dimensions();

	try {
		ismrmrd_dataset_->appendImageHeader(*m1->getObjectPtr(),"image.head");

		hoNDArray<std::complex<float> > *buffer_ = new hoNDArray<std::complex<float> >(dim, m2->getObjectPtr()->get_data_ptr(),false);
		
		std::vector<unsigned int> dims(dim.size());

		size_t i;
		for (i = 0; i< dim.size(); i++) {
			dims[i] = dim[i];
		}

		if (ismrmrd_dataset_->appendArray(dims, m2->getObjectPtr()->get_data_ptr(), "image_0.img") < 0) {
			GERROR("Failed to write image data\n");
			return GADGET_FAIL;
		}
	} catch (...) {
		GERROR("Error attempting to append images to HDF5 file\n");
		return GADGET_FAIL;
	}

	if (this->next()->putq(m1) == -1) {
		m1->release();
		ACE_ERROR_RETURN((LM_ERROR, ACE_TEXT("%p\n"), ACE_TEXT("ImageSaveHDFGadget::process, passing on data on to next gadget")), -1);
	}
	
	return GADGET_OK;
}

inline bool ImageSaveHDFGadgetCPLX::dirExists(const std::string &dirName_in)
{
	DWORD ftyp = GetFileAttributesA(dirName_in.c_str());

	if (ftyp == INVALID_FILE_ATTRIBUTES) {
		return false;		//something is wrong with your path!
	}

	if (ftyp & FILE_ATTRIBUTE_DIRECTORY) {
		return true;		// this is a directory!
	}

	return false;		// this is not a directory!
}

std::string ImageSaveHDFGadgetCPLX::get_date_time_string()
{
	time_t rawtime;
	struct tm *timeinfo;
	time (&rawtime);
	timeinfo = localtime (&rawtime);

	std::stringstream str;
	str << timeinfo->tm_year+1900
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_mon+1
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_mday
		<< "-"
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_hour
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_min
		<< std::setw(2) << std::setfill('0') << timeinfo->tm_sec;

	std::string ret = str.str();

	return ret;
}

GADGET_FACTORY_DECLARE(ImageSaveHDFGadgetCPLX)
