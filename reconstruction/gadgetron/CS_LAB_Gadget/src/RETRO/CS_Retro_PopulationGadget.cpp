#include "CS_Retro_PopulationGadget.h"

namespace Gadgetron{
// class constructor
CS_Retro_PopulationGadget::CS_Retro_PopulationGadget() {}
 
// class destructor - delete temporal buffer/memory
CS_Retro_PopulationGadget::~CS_Retro_PopulationGadget(){
	
}

// read flexible data header
int CS_Retro_PopulationGadget::process_config(ACE_Message_Block* mb)
{
	return GADGET_OK;
}

int CS_Retro_PopulationGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,GadgetContainerMessage< hoNDArray< float > >* m2, GadgetContainerMessage< hoNDArray <std::complex<float> > >* m3){
	int iNoGates = GlobalVar::instance()->iNoGates_;
	fTolerance_ = 2;
	iNoChannels_ = m3->getObjectPtr()->get_size(2);
	
	// get number of PE/PA lines
	vPE_ = GlobalVar::instance()->vPE_;
	vPA_ = GlobalVar::instance()->vPA_;
	
	// get navigator and convert to std::vector
	//hafNav_ = *m2->getObjectPtr();
	for (int iI = 0; iI < m2->getObjectPtr()->get_size(0); iI++)
		vNavInt_.push_back(m2->getObjectPtr()->at(iI));//hafNav_(iI));
	//GADGET_DEBUG2("global PE: %i, PA: %i\n", vPE_.size(), vPA_.size());
	// get unordered kspace data
	vtDims_unordered_ = *m3->getObjectPtr()->get_dimensions();
	hacfKSpace_unordered_.create(m3->getObjectPtr()->get_dimensions());
	memcpy(hacfKSpace_unordered_.get_data_ptr(), m3->getObjectPtr()->get_data_ptr(), sizeof(std::complex<float>)*m3->getObjectPtr()->get_number_of_elements());

	//-------------------------------------------------------------------------
	// discard first seconds of the acquisitions and wait for steady-state
	//-------------------------------------------------------------------------
	if (fDiscard()){
		GADGET_DEBUG1("Error occured in CS_Retro_Population::fDiscard(..) - process aborted\n");
		return GADGET_FAIL;
	}
	
	//-------------------------------------------------------------------------
	// get centroids
	//-------------------------------------------------------------------------
	if (fCalcCentroids(iNoGates)){
		GADGET_DEBUG1("Error occured in CS_Retro_Population::fCalcCentroids(..) - process aborted\n");
		return GADGET_FAIL;
	}
	else{
		for (int i = 0; i < vfCentroids_.size(); i++){
			GADGET_DEBUG2("Centroid %i: %f\n", i, vfCentroids_.at(i));
		}
	}

	//-------------------------------------------------------------------------
	// populate k-space: mode: closest, gates: 4
	//-------------------------------------------------------------------------
	if (fPopulatekSpace(iNoGates)){
		GADGET_DEBUG1("Error occured in CS_Retro_Population::fPopulatekSpace() - process aborted\n");
	}

	GADGET_DEBUG1("kSpace populated and ready to stream..\n");
	hacfKSpace_reordered_.print(std::cout);

	//-------------------------------------------------------------------------
	// make new GadgetContainer
	//-------------------------------------------------------------------------
	// create container
	GadgetContainerMessage<hoNDArray<std::complex<float>>>* tmp_m2 = new GadgetContainerMessage<hoNDArray<std::complex<float>>>();
	
	// concat
	m1->cont(tmp_m2);
	
	// create
	try{ tmp_m2->getObjectPtr()->create(hacfKSpace_reordered_.get_dimensions());}
	catch (std::runtime_error &err){
		GADGET_DEBUG_EXCEPTION(err, "CS_Retro: Unable to allocate new image array\m");
		m1->release();
		return -1;
	}
	
	memcpy(tmp_m2->getObjectPtr()->get_data_ptr(), hacfKSpace_reordered_.get_data_ptr(), sizeof(std::complex<float>)*hacfKSpace_reordered_.get_number_of_elements());
	
	// put on q
	if (this->next()->putq(m1) < 0) {
    	return GADGET_FAIL;
	}

	return GADGET_OK;
}

bool CS_Retro_PopulationGadget::fDiscard(){
	float fPreCrop = 5; // [s]
	int iStartIndex = std::floor(fPreCrop/(fTR_/1000));
	vNavInt_.erase(vNavInt_.begin(), vNavInt_.begin()+iStartIndex);
	vPA_.erase(vPA_.begin(), vPA_.begin() + iStartIndex);
	vPE_.erase(vPE_.begin(), vPE_.begin() + iStartIndex);	

	GADGET_DEBUG2("first seconds discarded - %i samples erased - TR: %f..\n", iStartIndex, fTR_);

	// new array size
	std::vector<size_t> vtDims_new = *hacfKSpace_unordered_.get_dimensions();
	vtDims_new.at(1) = vtDims_new.at(1) - iStartIndex;

	GADGET_DEBUG1("kSpace before deletion\n");
	hacfKSpace_unordered_.print(std::cout);

	// new array
	hoNDArray<std::complex<float>> hacfTmp(vtDims_new); hacfTmp.fill((0.0,0.0));
	for (int iR = 0; iR < hacfKSpace_unordered_.get_size(0); iR++){
		for (int iL = iStartIndex; iL <  hacfKSpace_unordered_.get_size(1); iL++){
			for (int iC = 0; iC < hacfKSpace_unordered_.get_size(2); iC++){
				hacfTmp(iR, iL-iStartIndex, iC) = hacfKSpace_unordered_(iR, iL, iC);
			}
		}
	}

	GADGET_DEBUG1("kSpace deleted:\n");
	hacfTmp.print(std::cout);
	hacfKSpace_unordered_ = hacfTmp;

	return GADGET_OK;
}

bool CS_Retro_PopulationGadget::fCalcCentroids(int iNoGates){
	// get centroids
	float fNavMin, fNavMax;
	switch(iGatingMode_){
		
		// percentile
		case 0:
			if (bMatlab_){
//				mexPrintf("get inhale/exhale borders by 10th and 90th percentile..\n");mexEvalString("drawnow;");
			}
			else{
				GADGET_DEBUG1("get inhale/exhale borders by 10th and 90th percentile..\n");
			}

			fNavMin = vNavInt_.at(std::min_element(vNavInt_.begin(), vNavInt_.end())-vNavInt_.begin());
			fNavMax = vNavInt_.at(std::max_element(vNavInt_.begin(), vNavInt_.end())-vNavInt_.begin());
			
			if (bMatlab_){
//				mexPrintf("navigator min: %.1f, max: %.1f\n", fNavMin, fNavMax);mexEvalString("drawnow;");
			}
			else{
				GADGET_DEBUG2("navigator min: %.1f, max: %.1f\n", fNavMin, fNavMax);
			}
			if (fNavMin == fNavMax){
				vfCentroids_.push_back(fNavMin);
			}
			else{
				// get histogram
				int iNumberBins = 256;
				std::vector<size_t> histogram_ = std::vector<size_t>(iNumberBins);
			
				// init 0
				for (size_t i = 0; i < iNumberBins; i++) histogram_.at(i) = 0;
				for (unsigned long int i = 0; i < vNavInt_.size(); i++){
					size_t bin = static_cast<size_t>(std::floor(vNavInt_.at(i)/((fNavMax-fNavMin)/iNumberBins)));
					if (bin >= iNumberBins)	bin = iNumberBins - 1;
					if (bin < 0) bin = 0;
					histogram_.at(bin)++;
				}
			
				// find 90th percentile
				long long cumsum = 0;
				size_t counter = 0;
				while (cumsum < (.90*vNavInt_.size())) cumsum += (long long)(histogram_.at(counter++));
				int   i90p = counter;
				float f90p = counter*((fNavMax-fNavMin)/iNumberBins);
			
				// find 10th percentile
				counter = 0;
				cumsum = 0;
				while (cumsum < (.10*vNavInt_.size())) cumsum += (long long)(histogram_[counter++]);
				int   i10p = counter;
				float f10p = counter*((fNavMax-fNavMin)/iNumberBins);
				if (bMatlab_){
//					mexPrintf("get equally spaced gate position - 10th: %.2f, 90th: %.2f, nPhases: %i\n", f10p, f90p, iNoGates);mexEvalString("drawnow;");
				}
				else{
					GADGET_DEBUG2("get equally spaced gate position - 10th: %.2f, 90th: %.2f, nPhases: %i\n", f10p, f90p, iNoGates);
				}

				// eqully spaced gate positions
				float fDistance = (f90p-f10p)/(iNoGates-1);
				for (long iI = 0; iI < iNoGates; iI++){
					vfCentroids_.push_back(f10p + iI*fDistance);
					if (bMatlab_){					
//						mexPrintf("Cent %f: %.5f\n", iI, vfCentroids_.at(iI));mexEvalString("drawnow;");
					}
				}								

				// get tolerance of the gate positions
				float fTolerance = std::abs(vfCentroids_.at(0)-vfCentroids_.at(1))*fTolerance_/2;

				// fill tolerance vector
				for (int i = 0; i < iNoGates; i++) vTolerance_.push_back(fTolerance);
			}
			break;

		// k-means
		case 1:
			if (bMatlab_){
//				mexPrintf("reorder_kSpace: k-means gating is not implemented in this version!\n");mexEvalString("drawnow;");
			}
			else{
				GADGET_DEBUG1("reorder_kSpace: k-means gating is not implemented in this version!\n");
			}
			GADGET_FAIL;
			break;

		default:
			if (bMatlab_){
//				mexPrintf("reorder_kSpace: no gating mode specified!\n");mexEvalString("drawnow;");
			}
			else{
				GADGET_DEBUG1("reorder_kSpace: no gating mode specified!\n");
			}
			return GADGET_FAIL;
			break;
	}

	return GADGET_OK;
}

bool CS_Retro_PopulationGadget::fPopulatekSpace(int iNoGates){
	if (bMatlab_){
//		mexPrintf("--- populate k-space ---");mexEvalString("drawnow;");
	}
	else{
		GADGET_DEBUG1("--- populate k-space ---");
	}	

	// drecks mdh
	if (vPE_.size()>hacfKSpace_unordered_.get_size(1)){
		vPE_.pop_back();
		vPA_.pop_back();
	}
	
	float fDist = std::sqrt(0.065)/2;
	
	// distinguish population mode
	switch(iPopulationMode_){
		
		// closest
		case 0:
			if (bMatlab_){
//				mexPrintf("closest mode..\n");mexEvalString("drawnow;");
			}
			else{
				GADGET_DEBUG1("closest mode..\n");
			}

			// initialize output k-space array
			hacfKSpace_reordered_.create(dimensionsIn_.at(0)*2, dimensionsIn_.at(1), dimensionsIn_.at(2), iNoGates, iNoChannels_);
			if (bMatlab_){
//				mexPrintf("global PE: %i, PA: %i\n", vPE_.size(), vPA_.size());mexEvalString("drawnow;");
			}
			else{
				GADGET_DEBUG2("global PE: %i, PA: %i\n", vPE_.size(), vPA_.size());			
			}
			
			// loop over phases/gates
			for (int iPh = 0; iPh <  iNoGates; iPh++){

				// get weights
				std::vector<float> vWeights(vNavInt_.size());
				for (long i = 0; i < vWeights.size(); i++)
					vWeights.at(i) = abs(vNavInt_.at(i) - vfCentroids_.at(iPh));
								
				if (bMatlab_){
//					mexPrintf("weights calculated - phase: %i\n", iPh);mexEvalString("drawnow;");
				}
				else{
					GADGET_DEBUG2("weights calculated - phase: %i\n", iPh);
				}

				// loop over lines
				for (int iLine = 0; iLine < dimensionsIn_.at(1); iLine++){
						
					// loop over partitions
					for (int iPar = 0; iPar < dimensionsIn_.at(2); iPar++){
						
						// check if iLine was acquired and push Indices on vector
						std::vector<long> lIndices;
						for (long  i = 0; i < vPE_.size(); i++)
							if (vPE_.at(i) == iLine)
								lIndices.push_back(i);
						
						// check iPar of the found acquisitions
						std::vector<long> lIndices2;
						for (long n = 0; n < lIndices.size(); n++)
							if (vPA_.at(lIndices.at(n)) == iPar)
								lIndices2.push_back(lIndices.at(n));
						
						// if no index is in the vector --> continue
						if (lIndices2.size() > 0){				
							// get weights (if multiple lines were found)
							std::vector<float> vThisDist; vThisDist.clear();
							for (int i = 0; i < lIndices2.size(); i ++) vThisDist.push_back(vWeights.at(lIndices2.at(i)));
							
							// min distance
							int iIndexMinDist = lIndices2.at(std::min_element(vThisDist.begin(), vThisDist.end())-vThisDist.begin());
							float fValMinDist = vWeights.at(iIndexMinDist);
							
							if (fValMinDist > vTolerance_.at(iPh))// && (abs((float)iLine - (float)iEchoLine_)) > fDist*dimensionsIn_[1] || (abs((float)iPar - (float)iEchoPartition_)) > fDist*dimensionsIn_[2])
								continue;
							
							// save acquisition into k-space - loop over channels
							for (int c = 0; c < iNoChannels_; c++){
								size_t tOffset_ordered = iLine*dimensionsIn_.at(0)*2+iPar*dimensionsIn_.at(1)*dimensionsIn_.at(0)*2+iPh*dimensionsIn_.at(2)*dimensionsIn_.at(1)*dimensionsIn_.at(0)*2+c*iNoGates*dimensionsIn_.at(2)*dimensionsIn_.at(1)*dimensionsIn_.at(0)*2;
								size_t tOffset_unordered = c*vtDims_unordered_.at(0)*vtDims_unordered_.at(1) + iIndexMinDist*vtDims_unordered_.at(0);
								memcpy(hacfKSpace_reordered_.get_data_ptr() + tOffset_ordered, hacfKSpace_unordered_.get_data_ptr() + tOffset_unordered, sizeof(std::complex<float>)*dimensionsIn_.at(0)*2);
							}
						}
					}
				}
				if (bMatlab_){
//					mexPrintf("kspace populated - phase: %i\n", iPh);mexEvalString("drawnow;");
				}
				else{
					GADGET_DEBUG2("kspace populated - phase: %i\n", iPh);
				}
			}

			break;
		 
		// average
		case 1:
			GADGET_DEBUG1("reorder_kSpace: population mode 'average' not implemented in this version\n");
			return GADGET_FAIL;
			break;

		// collect
		case 2:
			GADGET_DEBUG1("reorder_kSpace: population mode 'collect' not implemented in this version\n");
			return GADGET_FAIL;
			break;

		default:
			GADGET_DEBUG1("reorder_kSpace: no population mode specified!\n");
			return GADGET_FAIL;
			break;
	}

	return GADGET_OK;
}

GADGET_FACTORY_DECLARE(CS_Retro_PopulationGadget)
}