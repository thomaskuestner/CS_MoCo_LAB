#include "CS_MRIImageToAcquisitionGadget.h"

using namespace Gadgetron;

// convert the image data set to individual acquisitions
int CS_MRIImageToAcquisitionGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2)
{
	// get dimensions of incoming data (x,y,z,t,c)
	vDims_ = *m2->getObjectPtr()->get_dimensions();

	// ---------------------------------------------------------------------------------
	// ------------------------------ conversion ---------------------------------------
	// ---------------------------------------------------------------------------------
	for (size_t iPhase = 0; iPhase < vDims_[3]; iPhase++) {
		for (size_t iPartition = 0; iPartition < vDims_[2]; iPartition++) {
			for (size_t iLine = 0; iLine < vDims_[1]; iLine++) {

				// new AcquisitionHeader
				GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_hdr_m1 = new GadgetContainerMessage<ISMRMRD::AcquisitionHeader>();

				// init acquisition header
				memset(GC_acq_hdr_m1->getObjectPtr(), 0, sizeof(ISMRMRD::AcquisitionHeader));

				// copy data from global stored header
				fCopyAcqHeader(GC_acq_hdr_m1, GlobalVar::instance()->AcqVec_.at(0));

				// correct some header information
				fCorrectHeader(GC_acq_hdr_m1, iLine, iPartition, iPhase);

				// create empty 1D array (COL*CHA)
				GadgetContainerMessage< hoNDArray< std::complex<float> > >* hacfTmp = new GadgetContainerMessage<hoNDArray< std::complex<float> > >();

				try {
					hacfTmp->getObjectPtr()->create(vDims_[0]*m1->getObjectPtr()->channels);
				} catch (std::runtime_error &err) {
					GEXCEPTION(err, "Unable to allocate new image array\n");

					hacfTmp->release();
					return -1;
				}

				// fill 1D array
				std::complex<float> *cfPtrOut = hacfTmp->getObjectPtr()->get_data_ptr(), *cfPtrIn = m2->getObjectPtr()->get_data_ptr();
				for (int iC = 0; iC < m1->getObjectPtr()->channels; iC++) {
					int iOffsetOut	= iC*vDims_[0];
					int iOffsetIn	= iLine*vDims_[0] + iPartition*vDims_[0]*vDims_[1] + iPhase*vDims_[0]*vDims_[1]*vDims_[2] + iC*vDims_[0]*vDims_[1]*vDims_[2]*vDims_[3];
					memcpy(cfPtrOut + iOffsetOut, cfPtrIn + iOffsetIn, sizeof(std::complex<float>)*vDims_[0]*m1->getObjectPtr()->channels);
				}

				// concat header to data
				GC_acq_hdr_m1->cont(hacfTmp);

				// concatenate
				if (this->next()->putq(GC_acq_hdr_m1) < 0) {
					return GADGET_FAIL;
				}
			}
		}
	}

	// clear AcquisitionHeader vector
	GlobalVar::instance()->AcqVec_.clear();

	return GADGET_OK;
}

// correct header information for filled full kSpace data
int CS_MRIImageToAcquisitionGadget::fCorrectHeader(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *GC_acq_hdr_m1, size_t iLine, size_t iPartition, size_t iPhase)
{
	// fill temporal acquisition header with values from first scan
	fCopyAcqHeader(GC_acq_hdr_m1, GlobalVar::instance()->AcqVec_.at(0));

	// set loop counter
	GC_acq_hdr_m1->getObjectPtr()->idx.kspace_encode_step_1 = iLine;
	GC_acq_hdr_m1->getObjectPtr()->idx.kspace_encode_step_2 = iPartition;
	GC_acq_hdr_m1->getObjectPtr()->idx.phase				= iPhase;
	
	// correct scan counter
	GC_acq_hdr_m1->getObjectPtr()->scan_counter = iLine + iPartition*vDims_[1] + iPhase*vDims_[1]*vDims_[2];

	// clear false flags
	GC_acq_hdr_m1->getObjectPtr()->flags = 0;

	// ---------------------------------------------------------------------------------
	// -------------------------- set several flags ------------------------------------
	// ---------------------------------------------------------------------------------
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	// first scan in partition encoding
	if (iPartition == 0) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2-1);
	}

	// last scan in partition encoding
	if (iPartition == vDims_[2]-1) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2-1);
	}
	
	// first scan in phase encoding
	if (iLine == 0) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1-1);
	}

	// last scan in phase encoding
	if (iLine == vDims_[1]-1) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1-1);
	}

	// last scan in measurement
	if (iLine == vDims_[1]-1 && iPartition == vDims_[2]-1 && iPhase == vDims_[3]-1) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ISMRMRD_ACQ_LAST_IN_MEASUREMENT-1);
	}
#else
	if (iPartition == 0) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_FIRST_IN_ENCODE_STEP2-1);
	}

	// last scan in partition encoding
	if (iPartition == vDims_[2]-1) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_LAST_IN_ENCODE_STEP2-1);
	}
	
	// first scan in phase encoding
	if (iLine == 0) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_FIRST_IN_ENCODE_STEP1-1);
	}

	// last scan in phase encoding
	if (iLine == vDims_[1]-1) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_LAST_IN_ENCODE_STEP1-1);
	}

	// last scan in measurement
	if (iLine == vDims_[1]-1 && iPartition == vDims_[2]-1 && iPhase == vDims_[3]-1) {
		GC_acq_hdr_m1->getObjectPtr()->flags = GC_acq_hdr_m1->getObjectPtr()->flags | 1<<(ISMRMRD::ACQ_LAST_IN_MEASUREMENT-1);
	}
#endif

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_MRIImageToAcquisitionGadget)
