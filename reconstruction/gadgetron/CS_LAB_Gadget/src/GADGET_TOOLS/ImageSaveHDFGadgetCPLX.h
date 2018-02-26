// from IsmrmrdDumpGadget
// from http://stackoverflow.com/questions/8233842/how-to-check-if-directory-exist-using-c-and-winapi

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include <complex>
#include <ismrmrd.h>
#include <fstream>

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	#include <dataset.h>
#else
	#include <ismrmrd_hdf5.h>
#endif

#ifdef __WIN32__
	#include <Shlwapi.h>
#endif

#include <CS_LAB_export.h>

namespace Gadgetron
{
	class EXPORTCSLAB ImageSaveHDFGadgetCPLX : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<std::complex<float> > >
	{
	public:
		ImageSaveHDFGadgetCPLX()
		{
		}

		~ImageSaveHDFGadgetCPLX()
		{
		}

		GADGET_DECLARE(ImageSaveHDFGadgetCPLX);

	protected:
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);
		int process_config(ACE_Message_Block *mb);

		bool dirExists(const std::string &dirName_in);

		// make string with date and time
		std::string get_date_time_string();

		// file prefix (common: Image / Acquisition)
		std::string file_prefix_;
		std::string file_path_;
		std::string ismrmrd_file_name_;

#ifdef __GADGETRON_VERSION_HIGHER_3_6__
		boost::shared_ptr<ISMRMRD::ISMRMRD_Dataset> ismrmrd_dataset_;
#else
		boost::shared_ptr<ISMRMRD::IsmrmrdDataset>	ismrmrd_dataset_;
#endif

		bool append_timestamp_;
	};
}
