// from IsmrmrdDumpGadget
// from http://stackoverflow.com/questions/8233842/how-to-check-if-directory-exist-using-c-and-winapi

#pragma once

#include "Gadget.h"
#include "hoNDArray.h"
#include <complex>
#include <ismrmrd.h>
#include <fstream>
#include <ismrmrd_hdf5.h>
#include <Shlwapi.h>
#include <CS_LAB_export.h>

namespace Gadgetron
{
	class EXPORTCSLAB ImageSaveHDFGadgetFLOAT : public Gadget2<ISMRMRD::ImageHeader, hoNDArray<float> >
	{
	public:
		ImageSaveHDFGadgetFLOAT()
		{
		}

		~ImageSaveHDFGadgetFLOAT()
		{
		}

		GADGET_DECLARE(ImageSaveHDFGadgetFLOAT);

	protected:
		int process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,GadgetContainerMessage<hoNDArray<float> > *m2);

		// get filename, ..
		int process_config(ACE_Message_Block* mb);

		bool dirExists(const std::string& dirName_in);

		// make string with date and time
		std::string get_date_time_string();

		// file prefix (common: Image / Acquisition)
		std::string file_prefix_;
		std::string file_path_;
		std::string ismrmrd_file_name_;
		boost::shared_ptr<ISMRMRD::IsmrmrdDataset>  ismrmrd_dataset_;
		bool append_timestamp_;
	};
}
