#include <string>
#include <fstream>

#include <ismrmrd.h>

#include "debug_functions.h"

using namespace Gadgetron;

bool Gadgetron::all_elements_zero(const hoNDArray<std::complex<float> > &array)
{
	for (size_t i = 0; i < array.get_number_of_elements(); i++)
		if (abs(array.at(i)) != 0) {
			return false;
		}

	return true;
}

// print AcquisitionHeader (for debugging purposes only)
void Gadgetron::print_acquisition_header(int counter, const ISMRMRD::AcquisitionHeader &h)
{
	static ISMRMRD::AcquisitionHeader last_header;
	static bool header_set = false;

	std::cout << "Header " << counter << ": ";
	if (!header_set) {
		// print all
		std::cout << "version=" << h.version << ", ";
		std::cout << "flags=" << h.flags << ", ";
		std::cout << "measurement_uid=" << h.measurement_uid << ", ";
		std::cout << "scan_counter=" << h.scan_counter << ", ";
		std::cout << "acquisition_time_stamp=" << h.acquisition_time_stamp << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++)
			std::cout << "physiology_time_stamp[" << i << "]=" << h.physiology_time_stamp[i] << ", ";

		std::cout << "number_of_samples=" << h.number_of_samples << ", ";
		std::cout << "available_channels=" << h.available_channels << ", ";
		std::cout << "active_channels=" << h.active_channels << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_CHANNEL_MASKS; i++)
			std::cout << "channel_mask[" << i << "]=" << h.channel_mask[i] << ", ";

		std::cout << "discard_pre=" << h.discard_pre << ", ";
		std::cout << "discard_post=" << h.discard_post << ", ";
		std::cout << "center_sample=" << h.center_sample << ", ";
		std::cout << "encoding_space_ref=" << h.encoding_space_ref << ", ";
		std::cout << "trajectory_dimensions=" << h.trajectory_dimensions << ", ";
		std::cout << "sample_time_us=" << h.sample_time_us << ", ";

		for (int i = 0; i < 3; i++)
			std::cout << "position[" << i << "]=" << h.position[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "read_dir[" << i << "]=" << h.read_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "phase_dir[" << i << "]=" << h.phase_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "slice_dir[" << i << "]=" << h.slice_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			std::cout << "patient_table_position[" << i << "]=" << h.patient_table_position[i] << ", ";

		std::cout << "idx.kspace_encode_step_1=" << h.idx.kspace_encode_step_1 << ", ";
		std::cout << "idx.kspace_encode_step_2=" << h.idx.kspace_encode_step_2 << ", ";
		std::cout << "idx.average=" << h.idx.average << ", ";
		std::cout << "idx.slice=" << h.idx.slice << ", ";
		std::cout << "idx.contrast=" << h.idx.contrast << ", ";
		std::cout << "idx.phase=" << h.idx.phase << ", ";
		std::cout << "idx.repetition=" << h.idx.repetition << ", ";
		std::cout << "idx.set=" << h.idx.set << ", ";
		std::cout << "idx.segment=" << h.idx.segment << ", ";
		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			std::cout << "idx.user[" << i << "]=" << h.idx.user[i] << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			std::cout << "user_int[" << i << "]=" << h.user_int[i] << ", ";
		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++)
			std::cout << "user_float[" << i << "]=" << h.user_float[i] << ", ";
	} else {
		// print differences
		if (h.version != last_header.version)
			std::cout << "version=" << h.version << ", ";
		if (h.flags != last_header.flags)
			std::cout << "flags=" << h.flags << ", ";
		if (h.measurement_uid != last_header.measurement_uid)
			std::cout << "measurement_uid=" << h.measurement_uid << ", ";
		if (h.scan_counter != last_header.scan_counter)
			std::cout << "scan_counter=" << h.scan_counter << ", ";
		if (h.acquisition_time_stamp != last_header.acquisition_time_stamp)
			std::cout << "acquisition_time_stamp=" << h.acquisition_time_stamp << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_PHYS_STAMPS; i++)
			if (h.physiology_time_stamp[i] != last_header.physiology_time_stamp[i])
				std::cout << "physiology_time_stamp[" << i << "]=" << h.physiology_time_stamp[i] << ", ";

		if (h.number_of_samples != last_header.number_of_samples)
			std::cout << "number_of_samples=" << h.number_of_samples << ", ";
		if (h.available_channels != last_header.available_channels)
			std::cout << "available_channels=" << h.available_channels << ", ";
		if (h.active_channels != last_header.active_channels)
			std::cout << "active_channels=" << h.active_channels << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_CHANNEL_MASKS; i++)
			if (h.channel_mask[i] != last_header.channel_mask[i])
				std::cout << "channel_mask[" << i << "]=" << h.channel_mask[i] << ", ";

		if (h.discard_pre != last_header.discard_pre)
			std::cout << "discard_pre=" << h.discard_pre << ", ";
		if (h.discard_post != last_header.discard_post)
			std::cout << "discard_post=" << h.discard_post << ", ";
		if (h.center_sample != last_header.center_sample)
			std::cout << "center_sample=" << h.center_sample << ", ";
		if (h.encoding_space_ref != last_header.encoding_space_ref)
			std::cout << "encoding_space_ref=" << h.encoding_space_ref << ", ";
		if (h.trajectory_dimensions != last_header.trajectory_dimensions)
			std::cout << "trajectory_dimensions=" << h.trajectory_dimensions << ", ";
		if (h.sample_time_us != last_header.sample_time_us)
			std::cout << "sample_time_us=" << h.sample_time_us << ", ";

		for (int i = 0; i < 3; i++)
			if (h.position[i] != last_header.position[i])
				std::cout << "position[" << i << "]=" << h.position[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.read_dir[i] != last_header.read_dir[i])
				std::cout << "read_dir[" << i << "]=" << h.read_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.phase_dir[i] != last_header.phase_dir[i])
				std::cout << "phase_dir[" << i << "]=" << h.phase_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.slice_dir[i] != last_header.slice_dir[i])
				std::cout << "slice_dir[" << i << "]=" << h.slice_dir[i] << ", ";
		for (int i = 0; i < 3; i++)
			if (h.patient_table_position[i] != last_header.patient_table_position[i])
				std::cout << "patient_table_position[" << i << "]=" << h.patient_table_position[i] << ", ";

		if (h.idx.kspace_encode_step_1 != last_header.idx.kspace_encode_step_1)
			std::cout << "idx.kspace_encode_step_1=" << h.idx.kspace_encode_step_1 << ", ";
		if (h.idx.kspace_encode_step_2 != last_header.idx.kspace_encode_step_2)
			std::cout << "idx.kspace_encode_step_2=" << h.idx.kspace_encode_step_2 << ", ";
		if (h.idx.average != last_header.idx.average)
			std::cout << "idx.average=" << h.idx.average << ", ";
		if (h.idx.slice != last_header.idx.slice)
			std::cout << "idx.slice=" << h.idx.slice << ", ";
		if (h.idx.contrast != last_header.idx.contrast)
			std::cout << "idx.contrast=" << h.idx.contrast << ", ";
		if (h.idx.phase != last_header.idx.phase)
			std::cout << "idx.phase=" << h.idx.phase << ", ";
		if (h.idx.repetition != last_header.idx.repetition)
			std::cout << "idx.repetition=" << h.idx.repetition << ", ";
		if (h.idx.set != last_header.idx.set)
			std::cout << "idx.set=" << h.idx.set << ", ";
		if (h.idx.segment != last_header.idx.segment)
			std::cout << "idx.segment=" << h.idx.segment << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			if (h.idx.user[i] != last_header.idx.user[i])
				std::cout << "idx.user[" << i << "]=" << h.idx.user[i] << ", ";

		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_INTS; i++)
			if (h.user_int[i] != last_header.user_int[i])
				std::cout << "user_int[" << i << "]=" << h.user_int[i] << ", ";
		for (int i = 0; i < ISMRMRD::ISMRMRD_USER_FLOATS; i++)
			if (h.user_float[i] != last_header.user_float[i])
				std::cout << "user_float[" << i << "]=" << h.user_float[i] << ", ";
	}

	std::cout << std::endl;

	// set header
	last_header = h;
	header_set = true;
}

void Gadgetron::print_non_zero_elements(const hoNDArray<std::complex<float> > &array, const std::string &name, const int line)
{
	static int counter = 0;
	std::ofstream file_out;
	file_out.open("/tmp/"+name+"_"+std::to_string(line)+"_"+std::to_string(counter)+".txt");

	if (all_elements_zero(array)) {
		file_out << "All elements are zero!" << std::endl;
	} else {
		for (size_t i = 0; i < array.get_number_of_elements(); i++) {
			if (abs(array.at(i)) != 0) {
				file_out << i << ".: " << array.at(i) << std::endl;
			}
		}
	}

	file_out.close();
	counter++;
}

void Gadgetron::print_vector(const std::vector<float> &v)
{
	std::ofstream file_out;
	file_out.open("/tmp/vector_out.txt");

	for (std::vector<float>::const_iterator i = v.begin(); i != v.end(); ++i) {
		file_out << *i << std::endl;
	}

	file_out.close();
}

void Gadgetron::load_vector(std::vector<float> &v)
{
	v.clear();

	std::ifstream file_in("/tmp/vector_in.txt");

	float value;
	while (file_in >> value) {
		v.push_back(value);
	}
}

void Gadgetron::print_array(const hoNDArray<std::complex<float> > &array)
{
	std::ofstream file_head;
	file_head.open("/tmp/array_head.txt");

	array.print(file_head);

	file_head.close();

	std::ofstream file_out;
	file_out.open("/tmp/array_content.txt");

	for (size_t i = 0; i < array.get_number_of_elements(); i++) {
		file_out << array.at(i).real() << ", " << array.at(i).imag() << std::endl;
	}

	file_out.close();
}

void Gadgetron::load_array(hoNDArray<std::complex<float> > &array, const std::string &name, const bool is_boolean)
{
	if (is_boolean) {
		std::ifstream file_real(std::string("/opt/data/")+name+std::string(".txt"));

		bool value;
		size_t counter = 0;
		while (counter < array.get_number_of_elements() && file_real >> value) {
			if (value) {
				array.at(counter) = std::complex<float>(1, 0);
			} else {
				array.at(counter) = std::complex<float>(0, 0);
			}
			counter++;
		}
	} else {
		std::ifstream file_real(std::string("/opt/data/")+name+std::string("_real.txt"));
		std::ifstream file_imag(std::string("/opt/data/")+name+std::string("_imag.txt"));

		float real, imag;
		size_t counter = 0;
		while (counter < array.get_number_of_elements() && file_real >> real && file_imag >> imag) {
			array.at(counter) = std::complex<float>(real, imag);
			counter++;
		}
	}
}

void Gadgetron::load_array(hoNDArray<float> &array)
{
	std::ifstream file_real("/tmp/array_in.txt");

	float value;
	size_t counter = 0;
	while (counter < array.get_number_of_elements() && file_real >> value) {
		array.at(counter) = value;
		counter++;
	}
}
