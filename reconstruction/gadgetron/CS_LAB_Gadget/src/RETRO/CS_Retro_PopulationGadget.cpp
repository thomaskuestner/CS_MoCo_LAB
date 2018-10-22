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
	GlobalVar::instance()->iPopulationMode_			= PopulationMode.value();
	GlobalVar::instance()->cardiac_gating_mode_		= CardiacGatingMode.value();
	GlobalVar::instance()->respiratory_gating_mode_	= RespiratoryGatingMode.value();
	cardiac_tolerance_parameter_					= CardiacTolerance.value();
	respiratory_tolerance_parameter_				= RespiratoryTolerance.value();
	low_res_vs_										= LowResVS.value();
	omit_center_vs_									= OmitCenterVS.value();
#else
	GlobalVar::instance()->iPopulationMode_			= *(get_int_value("PopulationMode").get());
	GlobalVar::instance()->cardiac_gating_mode_		= *(get_int_value("CardiacGatingMode").get());
	GlobalVar::instance()->respiratory_gating_mode_	= *(get_int_value("RespiratoryGatingMode").get());
	cardiac_tolerance_parameter_					= *(get_float_value("CardiacTolerance").get());
	respiratory_tolerance_parameter_				= *(get_float_value("RespiratoryTolerance").get());
	low_res_vs_										= *(get_float_value("LowResVS").get());;
	omit_center_vs_									= *(get_float_value("OmitCenterVS").get());;
#endif

	return GADGET_OK;
}

int CS_Retro_PopulationGadget::process(GadgetContainerMessage<ISMRMRD::ImageHeader> *m1,GadgetContainerMessage<hoNDArray<float> > *m2, GadgetContainerMessage<hoNDArray<std::complex<float> > > *m3)
{
	iNoChannels_ = m3->getObjectPtr()->get_size(2);

	// get number of phases/gates
	unsigned int number_of_respiratory_phases = get_number_of_gates(m1->getObjectPtr()->user_int[0], 0);
	unsigned int number_of_cardiac_phases = get_number_of_gates(m1->getObjectPtr()->user_int[0], 1);

	// get navigator and convert to std::vector
	for (unsigned int iI = 0; iI < m2->getObjectPtr()->get_number_of_elements();) {
		navigator_resp_interpolated_.push_back(m2->getObjectPtr()->at(iI++));
		navigator_card_interpolated_.push_back(m2->getObjectPtr()->at(iI++));
	}

	// discard empty values from the end
	if (number_of_respiratory_phases == 1) {
		navigator_resp_interpolated_.clear();
	} else {
		discard_empty_elements_from_back(navigator_resp_interpolated_);
	}

	if (number_of_cardiac_phases == 1) {
		navigator_card_interpolated_.clear();
	} else {
		discard_empty_elements_from_back(navigator_card_interpolated_);
	}

	// correct number of phases if necessary
	if (navigator_card_interpolated_.size() == 0) {
		GWARN("No cardiac navigator signal found. Cardiac motion correction is disabled.\n");
		number_of_cardiac_phases = 1;
		set_number_of_gates(m1->getObjectPtr()->user_int[0], 1, 1);
	}

	if (navigator_resp_interpolated_.size() == 0) {
		GWARN("No respiratory navigator signal found. Respiratory motion correction is disabled.\n");
		number_of_respiratory_phases = 1;
		set_number_of_gates(m1->getObjectPtr()->user_int[0], 0, 1);
	}

	// get unordered kspace data
	hacfKSpace_unordered_.create(m3->getObjectPtr()->get_dimensions());
	memcpy(hacfKSpace_unordered_.get_data_ptr(), m3->getObjectPtr()->get_data_ptr(), sizeof(std::complex<float>)*m3->getObjectPtr()->get_number_of_elements());

	// free memory (also frees m3)
	m2->release();

	// initialize output k-space array (ReadOut x PhaseEncoding x PArtitionencoding x RespiratoryGates x CardiacGates x Channels)
	hacfKSpace_reordered_.create(hacfKSpace_unordered_.get_size(0), m1->getObjectPtr()->matrix_size[1], m1->getObjectPtr()->matrix_size[2], number_of_respiratory_phases, number_of_cardiac_phases, iNoChannels_);

	//-------------------------------------------------------------------------
	// discard first seconds of the acquisitions and wait for steady-state
	//-------------------------------------------------------------------------
	if (!fDiscard()) {
		GERROR("process aborted\n");
		return GADGET_FAIL;
	}

	//-------------------------------------------------------------------------
	// get respiratory centroids
	//-------------------------------------------------------------------------
	if (number_of_respiratory_phases > 1) {
		if (!get_respiratory_gates(number_of_respiratory_phases)) {
			GERROR("process aborted\n");
			return GADGET_FAIL;
		} else {
			for (size_t i = 0; i < respiratory_centroids_.size(); i++) {
				GDEBUG("Respiratory Centroid %i: %f\n", i, respiratory_centroids_.at(i));
			}
		}
	}

	//-------------------------------------------------------------------------
	// get cardiac centroids
	//-------------------------------------------------------------------------
	if (number_of_cardiac_phases > 1) {
		// f_s = 1000 -> ms raster
		//%fs = 1000; % ECGInt_ms on ms raster
		if (!get_cardiac_gates(number_of_cardiac_phases, 1000)) {
			GERROR("process aborted\n");
			return GADGET_FAIL;
		}
	}

	//-------------------------------------------------------------------------
	// populate k-space
	//-------------------------------------------------------------------------
	if (!fPopulatekSpace(number_of_cardiac_phases, number_of_respiratory_phases)) {
		GERROR("process aborted\n");
		return GADGET_FAIL;
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
	if (navigator_resp_interpolated_.size() >= static_cast<size_t>(iStartIndex)) {
		navigator_resp_interpolated_.erase(navigator_resp_interpolated_.begin(), navigator_resp_interpolated_.begin()+iStartIndex);
	} else {
		GWARN("more elements should be delete than there were actually in navigator_resp_interpolated_. Nothing is done!\n");
	}

	if (GlobalVar::instance()->vPA_.size() >= static_cast<size_t>(iStartIndex)) {
		GlobalVar::instance()->vPA_.erase(GlobalVar::instance()->vPA_.begin(), GlobalVar::instance()->vPA_.begin() + iStartIndex);
	} else {
		GWARN("more elements should be delete than there were actually in vPA_. Nothing is done!\n");
	}

	if (GlobalVar::instance()->vPE_.size() >= static_cast<size_t>(iStartIndex)) {
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
	hoNDArray<std::complex<float> > hacfTmp(vtDims_new);
	hacfTmp.fill(std::complex<float>(0.0, 0.0));
	#pragma omp parallel for
	for (size_t iR = 0; iR < hacfKSpace_unordered_.get_size(0); iR++) {
		for (size_t iL = static_cast<size_t>(iStartIndex); iL < hacfKSpace_unordered_.get_size(1); iL++) {
			for (size_t iC = 0; iC < hacfKSpace_unordered_.get_size(2); iC++) {
				hacfTmp(iR, iL-iStartIndex, iC) = hacfKSpace_unordered_(iR, iL, iC);
			}
		}
	}

	GINFO("kSpace deleted:\n");
	hacfTmp.print(std::cout);

	hacfKSpace_unordered_ = hacfTmp;

	return true;
}

bool CS_Retro_PopulationGadget::get_cardiac_gates(const unsigned int cardiac_gate_count, const float f_s)
{
	std::vector<size_t> x_pos;
	cardiac_gates_.clear();

	switch (GlobalVar::instance()->cardiac_gating_mode_) {
	// PanTompkin
	case 0:
		//%%Remove Outlier Peaks
		//%dEcg = fRemoveHighPeaks( dEcg, 'global',lDebug);
		//%[peak_r, loc_r, hfig_peak] = rpeak_detect(dEcg,fs,lDebug);
		remove_high_peaks(navigator_card_interpolated_);
		search_peaks(navigator_card_interpolated_, x_pos);

		break;

	// Wavelet
	case 1:
		GERROR("reorder_kSpace: Wavelet cardiac gating is not implemented in this version!\n");

		return false;
		break;

	// WaveletPanTomp
	case 2:
		GERROR("reorder_kSpace: WaveletPanTompkin cardiac gating is not implemented in this version!\n");

		return false;
		break;

	default:
		GERROR("reorder_kSpace: no cardiac gating mode specified!\n");

		return false;
		break;
	}

	//%% find RR-intervalls, HR & Ratio
	//%rr_int = diff(loc_r);
	std::vector<int> rr;// = arma::conv_to<std::vector<int> >::from(arma::diff(arma::Col<size_t>(x_pos)));

	//%% prevent a wrong RR interval detection
	//%lInd = rr_int >= 60/300 * fs; % hr is tob o large
	//%% rr_int = rr_int(lInd);
	//%lInd = find(lInd == 0)+1;
	//%loc_r(lInd) = [];
	//%peak_r(lInd) = [];
	//%rr_int = diff(loc_r);
	int deleted_elements = 0;
	for (size_t i = 0; i < rr.size(); i++) {
		if (rr.at(i) < 60.0/300.0 * f_s) {
			x_pos.erase(x_pos.begin() + i - deleted_elements + 1);
			deleted_elements++;
		}
	}
	rr = arma::conv_to<std::vector<int> >::from(arma::diff(arma::Col<size_t>(x_pos)));

	//%rr_mean = round(median(rr_int));
	int rr_mean = static_cast<int>(std::round(arma::median(arma::Col<int>(rr))));

	//%%Heartrate (1 / rr_int)
	//%hr = round(60./(rr_int./fs));
	std::vector<int> hr;
	for (size_t i = 0; i < rr.size(); i++) {
		hr.push_back(static_cast<int>(std::round(60.0/(rr.at(i)/f_s))));
	}

	switch (GlobalVar::instance()->cardiac_gating_mode_) {
	// linear
	case 0:
		//%dCardiacPhases = zeros(length(dEcg),1);
		//%for i=1:length(loc_r)-1
		//%	x = ceil(linspace(1,rr_int(i),iNcPhases+1));
		//%	y = [1:1:iNcPhases iNcPhases];
		//%	xy_sys = [1:rr_int(i)];
		//%	dCardiacPhases(loc_r(i):loc_r(i+1)-1) = floor(interp1(x,y,xy_sys,'linear'));
		//%end

		// fill up with zeros until the first interpolation
		for (size_t i = 1; i < x_pos.at(0); i++) {
			cardiac_gates_.push_back(0);
		}

		for (size_t i = 0; i < x_pos.size()-1; i++) {
			arma::vec x = arma::linspace(1, rr.at(i), cardiac_gate_count+1);
			for (size_t j = 0; j < x.n_elem; j++) {
				x.at(j) = std::ceil(x.at(j));
			}

			arma::vec y(cardiac_gate_count+1);
			for (unsigned int j = 0; j < cardiac_gate_count; j++) {
				y.at(j) = j+1;
			}
			y.at(y.n_elem-1) = cardiac_gate_count;

			arma::vec xy_sys(rr.at(i));
			for (int j = 0; j < rr.at(i); j++) {
				xy_sys.at(j) = j+1;
			}

			arma::vec interpolation;
			arma::interp1(x, y, xy_sys, interpolation, "*linear");

			for (size_t j = 0; j < interpolation.size(); j++) {
				cardiac_gates_.push_back(interpolation.at(j)-1);
			}
		}

		// fill up with zeros
		while (cardiac_gates_.size() < navigator_card_interpolated_.size()) {
			cardiac_gates_.push_back(0);
		}

		break;

	// diastole
	case 1:
		//%%specifiy systole/RR-ratio depending on heartrate
		//%rr_ratio = [linspace(0.33,0.33,60) linspace(0.33,0.5,61) linspace(0.5,0.5,120), linspace(0.5,0.5,max(hr)-241)];
		break;


	default:
		GERROR("reorder_kSpace: no cardiac gating mode specified!\n");

		return false;
		break;
	}

	if (cardiac_tolerance_parameter_ > 0.0) {
		// TODO: Make implementation
	}

	// TODO: implement rejection

	// TODO: implement crop

	return true;
}

bool CS_Retro_PopulationGadget::get_respiratory_gates(const unsigned int respiratory_gate_count)
{
	// get centroids
	float fNavMin, fNavMax;
	fNavMin = fNavMax = 0;

	switch (GlobalVar::instance()->respiratory_gating_mode_) {
	// percentile
	case 0:
		GINFO("get inhale/exhale borders by 10th and 90th percentile..\n");

		if (navigator_resp_interpolated_.size() > 0) {
			fNavMin = navigator_resp_interpolated_.at(std::min_element(navigator_resp_interpolated_.begin(), navigator_resp_interpolated_.end())-navigator_resp_interpolated_.begin());
			fNavMax = navigator_resp_interpolated_.at(std::max_element(navigator_resp_interpolated_.begin(), navigator_resp_interpolated_.end())-navigator_resp_interpolated_.begin());
		}

		GDEBUG("navigator min: %.1f, max: %.1f\n", fNavMin, fNavMax);

		if (fNavMin == fNavMax) {
			respiratory_centroids_.push_back(fNavMin);
		} else {
			// get histogram
			std::vector<unsigned int> histogram = std::vector<unsigned int>(256);

			for (size_t i = 0; i < navigator_resp_interpolated_.size(); i++) {
				int bin = static_cast<int>(std::floor((navigator_resp_interpolated_.at(i)-fNavMin)/((fNavMax-fNavMin)/histogram.size())));

				if (bin >= static_cast<int>(histogram.size())) {
					GWARN("Found bin %d >= %d. Corrected to bin = %d. navigator_resp_interpolated_ is %f\n", bin, histogram.size(), histogram.size()-1, navigator_resp_interpolated_.at(i));
					bin = histogram.size() - 1;
				}

				if (bin < 0) {
					GWARN("Found bin %d < 0. Corrected to bin = 0.\n", bin);
					bin = 0;
				}

				histogram.at(bin)++;
			}

			//%iSum = sum(iHist);
			unsigned int sum = std::accumulate(histogram.begin(), histogram.end(), 0, std::plus<unsigned int>());

			//%iLowerLim = 1;
			//%while sum(iHist(iLowerLim:end))/iSum > 0.9
			//%	iLowerLim = iLowerLim + 1;
			//%end
			//%iUpperLim = 256;
			//%while sum(iHist(1:iUpperLim))/iSum > 0.9
			//%	iUpperLim = iUpperLim - 1;
			//%end
			unsigned int lower_limit = 0;
			while (static_cast<float>(std::accumulate(histogram.begin()+lower_limit, histogram.end(), 0, std::plus<unsigned int>()))/static_cast<float>(sum) > 0.9) {
				lower_limit++;
			}

			unsigned int upper_limit = histogram.size()-1;
			while (static_cast<float>(std::accumulate(histogram.begin(), histogram.begin()+upper_limit, 0, std::plus<unsigned int>()))/static_cast<float>(sum) > 0.9) {
				upper_limit--;
			}

			//%dMin = dX(iLowerLim);
			//%dMax = dX(iUpperLim);
			const float stepsize = (fNavMax - fNavMin)/static_cast<float>(histogram.size());
			const float min = fNavMin + stepsize * lower_limit;
			const float max = fNavMin + stepsize * upper_limit;

			GINFO("get equally spaced gate position - 10th: %.2f, 90th: %.2f, respiratory phases: %i\n", min, max, respiratory_gate_count);

			//%dGatePos = linspace(dMax, dMin, iNPhasesOld)';
			respiratory_centroids_ = arma::conv_to<std::vector<float> >::from(arma::linspace(max, min, respiratory_gate_count));

			// get tolerance of the gate positions
			const float tolerance = std::abs(respiratory_centroids_.at(0)-respiratory_centroids_.at(1))*respiratory_tolerance_parameter_/2.0;

			// fill tolerance vector
			for (unsigned int i = 0; i < respiratory_gate_count; i++) {
				respiratory_tolerance_vector_.push_back(tolerance);
			}
		}
		break;

	// k-means
	case 1:
		GERROR("reorder_kSpace: k-means respiratory gating is not implemented in this version!\n");

		return false;
		break;

	default:
		GERROR("reorder_kSpace: no respiratory gating mode specified!\n");

		return false;
		break;
	}

	return true;
}

std::vector<float> CS_Retro_PopulationGadget::calculate_weights(const int population_mode, const unsigned int phase)
{
	std::vector<float> weights;

	if (respiratory_centroids_.size() <= phase) {
		return weights;
	} else {
		weights.resize(navigator_resp_interpolated_.size());
	}

	for (size_t i = 0; i < weights.size(); i++) {
		switch (population_mode) {
		// closest & average (weights are only needed for tolerance window)
		case 0:
		case 1:
			weights.at(i) = abs(navigator_resp_interpolated_.at(i) - respiratory_centroids_.at(phase));

			break;

		// gauss
		case 3:
			weights.at(i) = 1/(respiratory_tolerance_vector_.at(phase)*std::sqrt(2*M_PI)) * exp(-(std::pow(navigator_resp_interpolated_.at(i)-respiratory_centroids_.at(phase),2))/(2*(std::pow(respiratory_tolerance_vector_.at(phase),2))));
			break;

		default:
			throw runtime_error("Selected population_mode unknown!\n");
		}
	}

	return weights;
}

hoNDArray<std::complex<float> > CS_Retro_PopulationGadget::get_populated_data(const std::vector<size_t> &indices, const std::vector<float> &centroid_distances, int &cardiac_phase)
{
	hoNDArray<std::complex<float> > populated_data(hacfKSpace_unordered_.get_size(0), iNoChannels_);

	// Take shortcut by closest mode if it does not really matter (only 1 entry)
	int mode = GlobalVar::instance()->iPopulationMode_;
	if (centroid_distances.size() == 0 || indices.size() == 1) {
		mode = 0;
	}

	switch (mode) {
	// closest
	case 0: {
		// calculate min distance
		size_t index_min_this_dist;
		if (centroid_distances.size() > 0) {
			index_min_this_dist = std::min_element(centroid_distances.begin(), centroid_distances.end())-centroid_distances.begin();
		} else {
			index_min_this_dist = 0;
		}
		size_t index_min_dist = indices.at(index_min_this_dist);

		// copy data into return array
		size_t hacfKSpace_unordered__measurement_offset = index_min_dist * hacfKSpace_unordered_.get_size(0);
		#pragma omp parallel for
		for (size_t channel = 0; channel < populated_data.get_size(1); channel++) {
			size_t populated_channel_offset = channel * populated_data.get_size(0);
			size_t hacfKSpace_unordered__channel_offset = channel * hacfKSpace_unordered_.get_size(0) * hacfKSpace_unordered_.get_size(1);

			memcpy(populated_data.get_data_ptr() + populated_channel_offset, hacfKSpace_unordered_.get_data_ptr() + hacfKSpace_unordered__channel_offset + hacfKSpace_unordered__measurement_offset, sizeof(std::complex<float>)*populated_data.get_size(0));
		}

		// set cardiac phase
		if (cardiac_gates_.size() > index_min_dist) {
			cardiac_phase = cardiac_gates_.at(index_min_dist);
		} else {
			cardiac_phase = 0;
		}
	} break;

	// average
	case 1: {
		// pre-filling step with zeros may be omitted here because we set the data below (=) and do not just append (+=)

		// take all measurements into account
		for (size_t measurement = 0; measurement < indices.size(); measurement++) {
			size_t hacfKSpace_unordered__measurement_offset = indices.at(measurement) * hacfKSpace_unordered_.get_size(0);

			#pragma omp parallel for
			for (size_t channel = 0; channel < populated_data.get_size(1); channel++) {
				size_t populated_channel_offset = channel * populated_data.get_size(0);
				size_t hacfKSpace_unordered__channel_offset = channel * hacfKSpace_unordered_.get_size(0) * hacfKSpace_unordered_.get_size(1);

				#pragma omp parallel for
				for (size_t kx_pixel = 0; kx_pixel < populated_data.get_size(0); kx_pixel++) {
					// append weighted value to return array
					populated_data.at(kx_pixel + populated_channel_offset) = hacfKSpace_unordered_.at(kx_pixel + hacfKSpace_unordered__channel_offset + hacfKSpace_unordered__measurement_offset);
				}
			}
		}

		// calculate overall weight sum and divide by it (weighted addition must be 1 at end)
		float normalization_factor = indices.size();

		// only divide if possible
		if (normalization_factor != 0) {
			// divide all elements by weight sum
			populated_data /= normalization_factor;
		}

		// set cardiac phase
		// fill histogram
		size_t cardiac_phases[256] = {0};	// assume max 256 phases
		for (size_t measurement = 0; measurement < indices.size(); measurement++) {
			// increase count for present phase
			cardiac_phases[cardiac_gates_.at(indices.at(measurement))]++;
		}

		// choose maximum
		cardiac_phase = *std::max_element(cardiac_phases, cardiac_phases+ARRAYSIZE(cardiac_phases));
	} break;

	// gauss
	case 3: {
		// pre-fill return array with zeros (IMPORTANT STEP)
		populated_data.fill(std::complex<float>(0.0));

		// take all measurements into account
		for (size_t measurement = 0; measurement < indices.size(); measurement++) {
			size_t hacfKSpace_unordered__measurement_offset = indices.at(measurement) * hacfKSpace_unordered_.get_size(0);

			#pragma omp parallel for
			for (size_t channel = 0; channel < populated_data.get_size(1); channel++) {
				size_t populated_channel_offset = channel * populated_data.get_size(0);
				size_t hacfKSpace_unordered__channel_offset = channel * hacfKSpace_unordered_.get_size(0) * hacfKSpace_unordered_.get_size(1);

				#pragma omp parallel for
				for (size_t kx_pixel = 0; kx_pixel < populated_data.get_size(0); kx_pixel++) {
					// append weighted value to return array
					populated_data.at(kx_pixel + populated_channel_offset) += centroid_distances.at(measurement)*hacfKSpace_unordered_.at(kx_pixel + hacfKSpace_unordered__channel_offset + hacfKSpace_unordered__measurement_offset);
				}
			}
		}

		// calculate overall weight sum and divide by it (weighted addition must be 1 at end)
		float normalization_factor = std::accumulate(std::begin(centroid_distances), std::end(centroid_distances), 0.0, std::plus<float>());

		// only divide if possible
		if (normalization_factor != 0) {
			// divide all elements by weight sum
			populated_data /= normalization_factor;
		}

		// set cardiac phase
		// TODO: set phase
	} break;

	default:
		throw runtime_error("Selected population_mode unknown!\n");
	}

	return populated_data;
}

bool CS_Retro_PopulationGadget::fPopulatekSpace(const unsigned int cardiac_gate_count, const unsigned int respiratory_gate_count)
{
	GINFO("--- populate k-space ---\n");

	// check population mode
	switch (GlobalVar::instance()->iPopulationMode_) {
	case 0:
		GINFO("Using closest mode\n");
		break;
	case 1:
		GINFO("Using average mode\n");
		break;
	case 2:
		GERROR("reorder_kSpace: population mode 'collect' not implemented in this version\n");

		return false;
		break;
	case 3:
		GINFO("Using gauss mode\n");
		break;
	default:
		GERROR("reorder_kSpace: no population mode specified!\n");

		return false;
		break;
	}

	// drecks mdh
	if (GlobalVar::instance()->vPE_.size() > hacfKSpace_unordered_.get_size(1)) {
		GDEBUG("A vPE_/vPA_ measurement is removed!\n");
		GlobalVar::instance()->vPE_.pop_back();
		GlobalVar::instance()->vPA_.pop_back();
	}

	GDEBUG("global PE: %i, PA: %i\n", GlobalVar::instance()->vPE_.size(), GlobalVar::instance()->vPA_.size());

	// loop over phases/gates
	//#pragma omp parallel for (parallelizing line loop is more efficient and keeps output prints in order)
	for (unsigned int respiratory_phase = 0; respiratory_phase < respiratory_gate_count; respiratory_phase++) {
		// get weights
		std::vector<float> vWeights = calculate_weights(GlobalVar::instance()->iPopulationMode_, respiratory_phase);

		GINFO("weights calculated - phase: %i\n", respiratory_phase);

		// declare variable for later usage
		int cardiac_phase;

		// loop over lines
		#pragma omp parallel for
		for (size_t iLine = 0; iLine < hacfKSpace_reordered_.get_size(1); iLine++) {
			// loop over partitions
			#pragma omp parallel for
			for (size_t iPar = 0; iPar < hacfKSpace_reordered_.get_size(2); iPar++) {
				// check if iLine was acquired and push Indices on vector
				std::vector<size_t> lIndices;
				for (size_t i = 0; i < GlobalVar::instance()->vPE_.size(); i++) {
					if (GlobalVar::instance()->vPE_.at(i) == iLine) {
						lIndices.push_back(i);
					}
				}

				// check iPar of the found acquisitions
				std::vector<size_t> lIndices2;
				for (size_t n = 0; n < lIndices.size(); n++) {
					if (GlobalVar::instance()->vPA_.at(lIndices.at(n)) == iPar) {
						// only take measurement into account when tolerance is okay
						if (vWeights.size() == 0 || vWeights.at(lIndices.at(n)) < respiratory_tolerance_vector_.at(respiratory_phase)) {
							lIndices2.push_back(lIndices.at(n));
						}
					}
				}

				// if no index is in the vector --> continue
				if (lIndices2.size() > 0) {
					// get weights (if multiple lines were found)
					std::vector<float> vThisDist;

					if (vWeights.size() > 0) {
						for (size_t i = 0; i < lIndices2.size(); i++) {
							vThisDist.push_back(vWeights.at(lIndices2.at(i)));
						}
					}

					// populate the data
					hoNDArray<std::complex<float> > populated_data = get_populated_data(lIndices2, vThisDist, cardiac_phase);

					// correct cardiac phase if necessary
					if (cardiac_phase < 0) {
						GERROR("cardiac_phase (=%d) < 0. This should not happen!\n", cardiac_phase);
						cardiac_phase = 0;
					}

					// copy populated data into great reordered kspace
					#pragma omp parallel for
					for (int c = 0; c < iNoChannels_; c++) {
						size_t tOffset_reordered =
							iLine * hacfKSpace_reordered_.get_size(0)
							+ iPar * hacfKSpace_reordered_.get_size(1) * hacfKSpace_reordered_.get_size(0)
							+ respiratory_phase * hacfKSpace_reordered_.get_size(2) * hacfKSpace_reordered_.get_size(1) * hacfKSpace_reordered_.get_size(0)
							+ cardiac_phase * hacfKSpace_reordered_.get_size(3) * hacfKSpace_reordered_.get_size(2) * hacfKSpace_reordered_.get_size(1) * hacfKSpace_reordered_.get_size(0)
							+ c * hacfKSpace_reordered_.get_size(4) * hacfKSpace_reordered_.get_size(3) * hacfKSpace_reordered_.get_size(2) * hacfKSpace_reordered_.get_size(1) * hacfKSpace_reordered_.get_size(0);

						size_t offset_populated = c * populated_data.get_size(0);

						memcpy(hacfKSpace_reordered_.get_data_ptr() + tOffset_reordered, populated_data.get_data_ptr() + offset_populated, sizeof(std::complex<float>)*populated_data.get_size(0));
					}
				}
			}
		}

		GINFO("kspace populated - cardiac phase: %d, respiratory phase: %d\n", cardiac_phase, respiratory_phase);
	}

	return true;
}

template <typename T>
void CS_Retro_PopulationGadget::remove_high_peaks(std::vector<T> &signal)
{
	//%% Check if the ecg signal has a bias and remove it
	//%if (abs(mean(ecg)) > 1e-4)
	//%	ecg = ecg - mean(ecg);
	//%end
	remove_signal_bias(signal);

	//%globalAK(1): Multiplier of the standard deviation
	//%globalAK(2): Number of searching iterations
	//globalAK = [4.5 2];
	const float global_ak_1 = 4.5;
	const unsigned int global_ak_2 = 2;

	//%pIdx = zeros(length(ecg),1);
	std::vector<bool> p_idx;
	for (size_t i = 0; i < signal.size(); i++) {
		p_idx.push_back(false);
	}

	for (unsigned int i = 0; i < global_ak_2; i++) {
		//%globalTh = (globalAK(1)) *std( ecgRem );
		std::vector<T> global_th = signal;
		float std_dev = arma::stddev(arma::Col<T>(signal));
		std::transform(std::begin(global_th), std::end(global_th), std::begin(global_th), bind2nd(std::multiplies<T>(), global_ak_1*std_dev));

		//%pIdx = pIdx | (abs(ecgRem) > globalTh);
		//%ecgRem( pIdx ) = globalTh.*sign(ecg(pIdx));
		#pragma omp parallel for
		for (size_t i = 0; i < signal.size(); i++) {
			p_idx.at(i) = p_idx.at(i) | (std::abs(signal.at(i) > global_th.at(i)));
			if (p_idx.at(i)) {
				signal.at(i) = global_th.at(i) * sgn(signal.at(i));
			}
		}
	}

	//%% Remove offset
	//%ecgRem = ecgRem-mean(ecgRem);
	remove_signal_bias(signal);
}

template <typename T>
void CS_Retro_PopulationGadget::search_peaks(const std::vector<T> &signal, std::vector<size_t> &x_pos)
{
	//%xdiff = diff(x);
	std::vector<T> diff_signal = arma::conv_to<std::vector<T> >::from(arma::diff(arma::Col<T>(signal)));
	std::vector<T> y_pos;

	//%for i=1:length(xdiff)
	//%	if(xdiff(i) <= 0)
	//%		if(i ~= 1)
	//%			if((xdiff(i-1)) >= 0)
	//%				% max found
	//%				k = k + 1;
	//%				xpos(k) = i;
	//%				ypos(k) = x(i);
	//%				sumy = sumy + ypos(k);
	//%				if(k == 1)
	//%					sumx = xpos(k);
	//%				else
	//%					sumx = sumx + (xpos(k) - xpos(k-1));
	//%				end
	//%			end
	//%		end
	//%	end
	//%end
	for (size_t i = 0; i < diff_signal.size(); i++) {
		if (diff_signal.at(i) <= 0) {
			if (i != 0) {
				if (diff_signal.at(i-1) >= 0) {
					// max found
					x_pos.push_back(i);
					y_pos.push_back(signal.at(i));
				}
			}
		}
	}

	//%averagey = sumy / length(xpos);
	//%averagey = averagey / 1.5;
	float average_y = static_cast<float>(std::accumulate(std::begin(y_pos), std::end(y_pos), static_cast<T>(0), std::plus<T>())) / y_pos.size();
	average_y /= 1.5;

	//%averagex = sumx / length(xpos);
	// beware: strange sum calculation ;)
	float average_x = static_cast<float>(std::accumulate(std::begin(x_pos), std::end(x_pos), static_cast<T>(0), std::plus<T>())) / x_pos.size();
// 	float average_x = static_cast<float>(x_pos.at(x_pos.size()-1)) / x_pos.size();

	//%z = 0;
	//%for i=1:length(xpos)-1
	//%	if(ypos(i-z) < averagey && (xpos(i+1-z) - xpos(i-z)) < averagex)
	//%		% max found
	//%		xpos(i-z) = []; % get rid of outlayers that are not part of the cardiac motion
	//%		ypos(i-z) = [];
	//%		z = z + 1;
	//%	end
	//%end
	// NOTE: in C++ implementation, i is counted differently, so that no z is needed
	size_t i = 0;
	while (i < x_pos.size()-1) {
		if (y_pos.at(i) < average_y && (x_pos.at(i+1) - x_pos.at(i) < average_x)) {
			// max found
			x_pos.erase(x_pos.begin()+i);
			y_pos.erase(y_pos.begin()+i);
		} else {
			i++;
		}
	}
}

GADGET_FACTORY_DECLARE(CS_Retro_PopulationGadget)
