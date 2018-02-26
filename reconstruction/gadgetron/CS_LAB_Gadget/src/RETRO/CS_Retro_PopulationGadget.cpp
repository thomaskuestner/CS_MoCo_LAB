#include "CS_Retro_PopulationGadget.h"

using namespace Gadgetron;

// class constructor
CS_Retro_PopulationGadget::CS_Retro_PopulationGadget()
{
}

// class destructor - delete temporal buffer/memory
CS_Retro_PopulationGadget::~CS_Retro_PopulationGadget()
{
}

// read flexible data header
int CS_Retro_PopulationGadget::process_config(ACE_Message_Block *mb)
{
	// set properties
#ifdef __GADGETRON_VERSION_HIGHER_3_6__
	GlobalVar::instance()->iPopulationMode_		= PopulationMode.value();
	GlobalVar::instance()->iGatingMode_			= GatingMode.value();
#else
	GlobalVar::instance()->iPopulationMode_		= *(get_int_value("PopulationMode").get());
	GlobalVar::instance()->iGatingMode_			= *(get_int_value("GatingMode").get());
#endif

	return GADGET_OK;
}

int CS_Retro_PopulationGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,GadgetContainerMessage<hoNDArray<float> > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3)
{
	fTolerance_ = 2;
	iNoChannels_ = m3->getObjectPtr()->get_size(2);
	iNPhases_ = m1->getObjectPtr()->user_int[0];

	// get navigator and convert to std::vector
	//hafNav_ = *m2->getObjectPtr();
	for (int iI = 0; iI < m2->getObjectPtr()->get_size(0); iI++) {
		vNavInt_.push_back(m2->getObjectPtr()->at(iI));//hafNav_(iI));
	}

	// get unordered kspace data
	vtDims_unordered_ = *m3->getObjectPtr()->get_dimensions();

	// initialization of dimensionsIn_ (copy of AccumulatorGadget)
	// TODO: check implementation!
	dimensionsIn_[0] = m1->getObjectPtr()->matrix_size[0];
	dimensionsIn_[1] = m1->getObjectPtr()->matrix_size[1];
	dimensionsIn_[2] = m1->getObjectPtr()->matrix_size[2];

	hacfKSpace_unordered_.create(m3->getObjectPtr()->get_dimensions());
	memcpy(hacfKSpace_unordered_.get_data_ptr(), m3->getObjectPtr()->get_data_ptr(), sizeof(std::complex<float>)*m3->getObjectPtr()->get_number_of_elements());

	//-------------------------------------------------------------------------
	// discard first seconds of the acquisitions and wait for steady-state
	//-------------------------------------------------------------------------
	if (fDiscard()) {
		GERROR("process aborted\n");
		return GADGET_FAIL;
	}

	//-------------------------------------------------------------------------
	// get centroids
	//-------------------------------------------------------------------------
	if (fCalcCentroids(iNPhases_)) {
		GERROR("process aborted\n");
		return GADGET_FAIL;
	} else {
		for (int i = 0; i < vfCentroids_.size(); i++) {
			GDEBUG("Centroid %i: %f\n", i, vfCentroids_.at(i));
		}
	}

	//-------------------------------------------------------------------------
	// populate k-space: mode: closest, gates: 4
	//-------------------------------------------------------------------------
	if (fPopulatekSpace(iNPhases_)) {
		GERROR("process aborted\n");
	}

	GINFO("kSpace populated and ready to stream..\n");
	hacfKSpace_reordered_.print(std::cout);

	//-------------------------------------------------------------------------
	// make new GadgetContainer
	//-------------------------------------------------------------------------
	// create container
	GadgetContainerMessage<hoNDArray<std::complex<float> > > *tmp_m2 = new GadgetContainerMessage<hoNDArray<std::complex<float> > >();

	// concat
	m1->cont(tmp_m2);

	// create
	try{
		tmp_m2->getObjectPtr()->create(hacfKSpace_reordered_.get_dimensions());
		memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), hacfKSpace_reordered_.get_data_ptr(), sizeof(std::complex<float>)*hacfKSpace_reordered_.get_number_of_elements());
	} catch (std::runtime_error &err) {
		GEXCEPTION(err, "Unable to allocate new image array\n");
		m1->release();
		return GADGET_FAIL;
	}

	// put on q
	if (this->next()->putq(m1) < 0) {
		return GADGET_FAIL;
	}

	return GADGET_OK;
}

bool CS_Retro_PopulationGadget::fDiscard()
{
	float fPreCrop = 5; // [s]
	int iStartIndex = std::floor(fPreCrop/(GlobalVar::instance()->fTR_/1000));

	// Set iStartIndex to 0 when it is negative
	if (iStartIndex < 0) {
		GWARN("iStartIndex=%d < 0. It is set to 0, maybe you want to check your code?\n", iStartIndex);
		iStartIndex = 0;
	}

	// only erase data when there is enough
	if (vNavInt_.size() >= iStartIndex) {
		vNavInt_.erase(vNavInt_.begin(), vNavInt_.begin()+iStartIndex);
	} else {
		GWARN("more elements should be delete than there were actually in vNavInt_. Nothing is done!\n");
	}

	if (GlobalVar::instance()->vPA_.size() >= iStartIndex) {
		GlobalVar::instance()->vPA_.erase(GlobalVar::instance()->vPA_.begin(), GlobalVar::instance()->vPA_.begin() + iStartIndex);
	} else {
		GWARN("more elements should be delete than there were actually in vPA_. Nothing is done!\n");
	}

	if (GlobalVar::instance()->vPE_.size() >= iStartIndex) {
		GlobalVar::instance()->vPE_.erase(GlobalVar::instance()->vPE_.begin(), GlobalVar::instance()->vPE_.begin() + iStartIndex);
	} else {
		GWARN("more elements should be delete than there were actually in vPE_. Nothing is done!\n");
	}

	GINFO("first seconds discarded - %i samples erased - TR: %f..\n", iStartIndex, GlobalVar::instance()->fTR_);

	// new array size
	std::vector<size_t> vtDims_new = *hacfKSpace_unordered_.get_dimensions();
	vtDims_new.at(1) = vtDims_new.at(1) - iStartIndex;

	GINFO("kSpace before deletion\n");
	hacfKSpace_unordered_.print(std::cout);

	// new array
	hoNDArray<std::complex<float>> hacfTmp(vtDims_new); hacfTmp.fill((0.0,0.0));
	for (int iR = 0; iR < hacfKSpace_unordered_.get_size(0); iR++) {
		for (int iL = iStartIndex; iL < hacfKSpace_unordered_.get_size(1); iL++) {
			for (int iC = 0; iC < hacfKSpace_unordered_.get_size(2); iC++) {
				hacfTmp(iR, iL-iStartIndex, iC) = hacfKSpace_unordered_(iR, iL, iC);
			}
		}
	}

	GINFO("kSpace deleted:\n");
	hacfTmp.print(std::cout);

	hacfKSpace_unordered_ = hacfTmp;

	return GADGET_OK;
}

bool CS_Retro_PopulationGadget::fCalcCentroids(int iNoGates)
{
	// get centroids
	float fNavMin, fNavMax;
	fNavMin = fNavMax = 0;

	switch (GlobalVar::instance()->iGatingMode_) {
	// percentile
	case 0:
		GINFO("get inhale/exhale borders by 10th and 90th percentile..\n");

		if (vNavInt_.size() > 0) {
			fNavMin = vNavInt_.at(std::min_element(vNavInt_.begin(), vNavInt_.end())-vNavInt_.begin());
			fNavMax = vNavInt_.at(std::max_element(vNavInt_.begin(), vNavInt_.end())-vNavInt_.begin());
		}

		GDEBUG("navigator min: %.1f, max: %.1f\n", fNavMin, fNavMax);

		if (fNavMin == fNavMax) {
			vfCentroids_.push_back(fNavMin);
		} else {
			// get histogram
			int iNumberBins = 256;
			std::vector<size_t> histogram = std::vector<size_t>(iNumberBins);

			// init 0
			for (size_t i = 0; i < iNumberBins; i++) {
				histogram.at(i) = 0;
			}

			for (unsigned long int i = 0; i < vNavInt_.size(); i++) {
				size_t bin = static_cast<size_t>(std::floor(vNavInt_.at(i)/((fNavMax-fNavMin)/iNumberBins)));

				if (bin >= iNumberBins) {
					bin = iNumberBins - 1;
				}

				if (bin < 0) {
					bin = 0;
				}

				histogram.at(bin)++;
			}

			// find 90th percentile
			long long cumsum = 0;
			size_t counter = 0;
			while (cumsum < (.90*vNavInt_.size())) {
				cumsum += static_cast<long long>(histogram.at(counter++));
			}

			int		i90p = counter;
			float	f90p = counter*((fNavMax-fNavMin)/iNumberBins);

			// find 10th percentile
			counter = 0;
			cumsum = 0;
			while (cumsum < (.10*vNavInt_.size())) {
				cumsum += static_cast<long long>(histogram[counter++]);
			}

			int		i10p = counter;
			float	f10p = counter*((fNavMax-fNavMin)/iNumberBins);

			GINFO("get equally spaced gate position - 10th: %.2f, 90th: %.2f, nPhases: %i\n", f10p, f90p, iNoGates);

			// eqully spaced gate positions
			float fDistance = (f90p-f10p)/(iNoGates-1);
			for (long iI = 0; iI < iNoGates; iI++) {
				vfCentroids_.push_back(f10p + iI*fDistance);
			}

			// get tolerance of the gate positions
			float fTolerance = std::abs(vfCentroids_.at(0)-vfCentroids_.at(1))*fTolerance_/2.0;

			// fill tolerance vector
			for (int i = 0; i < iNoGates; i++) {
				vTolerance_.push_back(fTolerance);
			}
		}
		break;

	// k-means
	case 1:
		GERROR("reorder_kSpace: k-means gating is not implemented in this version!\n");

		return GADGET_FAIL;
		break;

	default:
		GERROR("reorder_kSpace: no gating mode specified!\n");

		return GADGET_FAIL;
		break;
	}

	return GADGET_OK;
}

bool CS_Retro_PopulationGadget::fPopulatekSpace(int iNoGates)
{
	GINFO("--- populate k-space ---\n");

	// drecks mdh
	if (GlobalVar::instance()->vPE_.size() > hacfKSpace_unordered_.get_size(1)) {
		GlobalVar::instance()->vPE_.pop_back();
		GlobalVar::instance()->vPA_.pop_back();
	}

	float fDist = std::sqrt(0.065)/2;

	// distinguish population mode
	switch(GlobalVar::instance()->iPopulationMode_) {
	// closest
	case 0:
		GINFO("closest mode..\n");

		// initialize output k-space array
		hacfKSpace_reordered_.create(dimensionsIn_[0]*2, dimensionsIn_[1], dimensionsIn_[2], iNoGates, iNoChannels_);

		GDEBUG("global PE: %i, PA: %i\n", GlobalVar::instance()->vPE_.size(), GlobalVar::instance()->vPA_.size());

		// loop over phases/gates
		for (int iPh = 0; iPh < iNoGates; iPh++) {
			// get weights
			std::vector<float> vWeights(vNavInt_.size());
			for (long i = 0; i < vWeights.size(); i++) {
				vWeights.at(i) = abs(vNavInt_.at(i) - vfCentroids_.at(iPh));
			}

			GINFO("weights calculated - phase: %i\n", iPh);

			// loop over lines
			for (int iLine = 0; iLine < dimensionsIn_[1]; iLine++) {
				// loop over partitions
				for (int iPar = 0; iPar < dimensionsIn_[2]; iPar++) {
					// check if iLine was acquired and push Indices on vector
					std::vector<long> lIndices;
					for (long i = 0; i < GlobalVar::instance()->vPE_.size(); i++) {
						if (GlobalVar::instance()->vPE_.at(i) == iLine) {
							lIndices.push_back(i);
						}
					}

					// check iPar of the found acquisitions
					std::vector<long> lIndices2;
					for (long n = 0; n < lIndices.size(); n++) {
						if (GlobalVar::instance()->vPA_.at(lIndices.at(n)) == iPar) {
							lIndices2.push_back(lIndices.at(n));
						}
					}

					// if no index is in the vector --> continue
					if (lIndices2.size() > 0) {
						// get weights (if multiple lines were found)
						std::vector<float> vThisDist;
						vThisDist.clear();

						for (int i = 0; i < lIndices2.size(); i ++) {
							vThisDist.push_back(vWeights.at(lIndices2.at(i)));
						}

						// min distance
						int iIndexMinDist = lIndices2.at(std::min_element(vThisDist.begin(), vThisDist.end())-vThisDist.begin());
						float fValMinDist = vWeights.at(iIndexMinDist);

						if (fValMinDist > vTolerance_.at(iPh)) {// && (abs((float)iLine - (float)iEchoLine_)) > fDist*dimensionsIn_[1] || (abs((float)iPar - (float)iEchoPartition_)) > fDist*dimensionsIn_[2])
							continue;
						}

						// save acquisition into k-space - loop over channels
						size_t max_offset_reordered, max_offset_unordered;
						max_offset_reordered = hacfKSpace_reordered_.get_number_of_elements() - 1;
						max_offset_unordered = hacfKSpace_unordered_.get_number_of_elements() - 1;
						for (int c = 0; c < iNoChannels_; c++) {
							size_t tOffset_reordered = iLine*dimensionsIn_[0]*2+iPar*dimensionsIn_[1]*dimensionsIn_[0]*2+iPh*dimensionsIn_[2]*dimensionsIn_[1]*dimensionsIn_[0]*2+c*iNoGates*dimensionsIn_[2]*dimensionsIn_[1]*dimensionsIn_[0]*2;
							size_t tOffset_unordered = c*vtDims_unordered_.at(0)*vtDims_unordered_.at(1) + iIndexMinDist*vtDims_unordered_.at(0);

							// protection against unallowed memory access
							if (tOffset_reordered > max_offset_reordered) {
								GWARN("Ordered offset is larger than allowed! current=%d, max=%d. Data will be corrupted!\n", tOffset_reordered, max_offset_reordered);
								break;
							}

							if (tOffset_unordered > max_offset_unordered) {
								GWARN("Unordered offset is larger than allowed! current=%d, max=%d. Data will be corrupted!\n", tOffset_unordered, max_offset_unordered);
								break;
							}

							memcpy(hacfKSpace_reordered_.get_data_ptr() + tOffset_reordered, hacfKSpace_unordered_.get_data_ptr() + tOffset_unordered, sizeof(std::complex<float>)*dimensionsIn_[0]*2);
						}
					}
				}
			}

			GINFO("kspace populated - phase: %i\n", iPh);
		}

		break;

	// average
	case 1:
		GERROR("reorder_kSpace: population mode 'average' not implemented in this version\n");

		return GADGET_FAIL;
		break;

	// collect
	case 2:
		GERROR("reorder_kSpace: population mode 'collect' not implemented in this version\n");

		return GADGET_FAIL;
		break;

	// gauss
	case 3:
		hacfKSpace_reordered_.create(dimensionsIn_[0]*2, dimensionsIn_[1], dimensionsIn_[2], iNoGates, iNoChannels_);

		//% CARDIAC PHASES loop
		//for icPh = iNPhasesCardLoop
		// loop over phases/gates
		for (int iPh = 0; iPh < iNoGates; iPh++) {
			//initialize the weighting vector
			std::vector<float> vWeights(vNavInt_.size());

			//fill Weightvector by Gauss-function
			// x = vNavInt_; sigma = vTolerance_ µ = cfVentorids_
			double dWeightAccu = 0;
			for (long k=0; k < vNavInt_.size(); k++) {
				vWeights.at(k) = 1/(vTolerance_.at(iPh)*std::sqrt(2*M_PI)) * exp(-(std::pow(vNavInt_.at(k)-vfCentroids_.at(iPh),2))/(2*(std::pow(vTolerance_.at(iPh),2))));
				dWeightAccu = dWeightAccu + vWeights.at(k);
			}

			if (dWeightAccu == 0) { //don't divide by 0
				dWeightAccu = 1;
			}

			GINFO("weights calculated - phase: %i\n", iPh);

			//	                dKSpace = complex(zeros(iNLines, iBaseRes.*2, iNPartitions, length(iNPhasesLoop), 1, iNChannels,sPrecision), ...
			//	                      zeros(iNLines, iBaseRes.*2, iNPartitions, length(iNPhasesLoop), 1, iNChannels,sPrecision));

			//	            lCardiacMask = sum(dCardiacPhases == icPh,2) > 0;
			//
			//	            % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
			//	            % PHASES loop
			//	            for iPh = iNPhasesLoop
			//	                lRespMask = dRespPhases == iPh;
			//	                dWeight = normpdf( dNavInt, dGatePos(iPh) , max(dTolerance(iPh,:))/sigma )';

			//	                for iL = 1:iNLines
			//	                    lLineMask = iLine == (iL - 1);
			// loop over lines
			for (int iLine = 0; iLine < dimensionsIn_[1]; iLine++) {
				// loop over partitions
				for (int iPar = 0; iPar < dimensionsIn_[2]; iPar++) {
					std::vector<long> lIndices;
					for (long i = 0; i < GlobalVar::instance()->vPE_.size(); i++) {
						if (GlobalVar::instance()->vPE_.at(i) == iLine) {
							lIndices.push_back(i);
						}
					}

					std::vector<long> lIndices2;
					for (long n = 0; n < lIndices.size(); n++) {
						if (GlobalVar::instance()->vPA_.at(lIndices.at(n)) == iPar) {
							lIndices2.push_back(lIndices.at(n));
						}
					}

					// if no index is in the vector --> continue
					if (lIndices2.size() > 0) {
						// get weights (if multiple lines were found)
						std::vector<float> vThisDist;
						vThisDist.clear();
						for (int i = 0; i < lIndices2.size(); i ++) {
							vThisDist.push_back(vWeights.at(lIndices2.at(i)));
							int currentindex = lIndices2.at(i);
							// min distance
							//int iIndexMinDist = lIndices2.at(std::min_element(vThisDist.begin(), vThisDist.end())-vThisDist.begin());
							//float fValMinDist = vWeights.at(iIndexMinDist);

							// save acquisition into k-space - loop over channels
							size_t max_offset_reordered, max_offset_unordered;
							max_offset_reordered = hacfKSpace_reordered_.get_number_of_elements() - 1;
							max_offset_unordered = hacfKSpace_unordered_.get_number_of_elements() - 1;
							for (int c = 0; c < iNoChannels_; c++) {
								size_t tOffset_reordered = iLine*dimensionsIn_[0]*2+iPar*dimensionsIn_[1]*dimensionsIn_[0]*2+iPh*dimensionsIn_[2]*dimensionsIn_[1]*dimensionsIn_[0]*2+c*iNoGates*dimensionsIn_[2]*dimensionsIn_[1]*dimensionsIn_[0]*2;
								size_t tOffset_unordered = c*vtDims_unordered_.at(0)*vtDims_unordered_.at(1) + currentindex*vtDims_unordered_.at(0);

								// protection against unallowed memory access
								if (tOffset_reordered > max_offset_reordered) {
									GWARN("Ordered offset is larger than allowed! current=%d, max=%d. Data will be corrupted!\n", tOffset_reordered, max_offset_reordered);
									break;
								}

								if (tOffset_unordered > max_offset_unordered) {
									GWARN("Unordered offset is larger than allowed! current=%d, max=%d. Data will be corrupted!\n", tOffset_unordered, max_offset_unordered);
									break;
								}

								memcpy(hacfKSpace_reordered_.get_data_ptr() + tOffset_reordered, hacfKSpace_unordered_.get_data_ptr() + tOffset_unordered, sizeof(std::complex<float>)*dimensionsIn_[0]*2);

								for(int x = 0; x < (dimensionsIn_[0]*2); x++) {
									// /dWeightAccu added (see loop below)
									// TODO: Error check here!
									hacfKSpace_reordered_.at(tOffset_reordered + x) = hacfKSpace_reordered_.at(tOffset_reordered + x) * vWeights.at(lIndices2.at(i)) / static_cast<float>(dWeightAccu);
								}
							}
						}
					}
				}

				//dDataAccu = dDataAccu./(dWeightAccu);
				// brought out to outer loop due to performance reasons.
				// correct equivalent formulation would be:
				// 		hacfKSpace_reordered_.at(i) = hacfKSpace_reordered_.at(i) / std::pow(static_cast<float>(dWeightAccu), (dimensionsIn_.at(1)*dimensionsIn_.at(2)));
				// Guess there's an error in that algorithm!
				// TODO: Error check here!
// 					for(long i = 0; i < hacfKSpace_reordered_.get_number_of_elements(); i++) {
// 						hacfKSpace_reordered_.at(i) = hacfKSpace_reordered_.at(i) / static_cast<float>(dWeightAccu);	// cast dWeightAccu to float in order to match operator /
// 					}
			}

			GINFO("kspace populated - phase: %i\n", iPh);
		}
		break;
	default:
		GERROR("reorder_kSpace: no population mode specified!\n");

		return GADGET_FAIL;
		break;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_PopulationGadget)
