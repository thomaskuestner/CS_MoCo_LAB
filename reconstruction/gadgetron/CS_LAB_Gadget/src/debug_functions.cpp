#include <string>
#include <fstream>

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

void Gadgetron::load_array(hoNDArray<std::complex<float> > &array)
{
	std::ifstream file_real("/tmp/array_real.txt");
	std::ifstream file_imag("/tmp/array_imag.txt");

	float real, imag;
	size_t counter = 0;
	while (counter < array.get_number_of_elements() && file_real >> real && file_imag >> imag) {
		array.at(counter) = std::complex<float>(real, imag);
		counter++;
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
