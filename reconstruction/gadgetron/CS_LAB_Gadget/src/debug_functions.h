#ifndef DEBUG_FUNCTIONS_H
#define DEBUG_FUNCTIONS_H

#include <string>
#include "hoNDArray.h"

#define print_hacfarray(name) print_non_zero_elements(name, "name", __LINE__)

namespace Gadgetron {
	bool all_elements_zero(const hoNDArray<std::complex<float> > &array);
	void print_non_zero_elements(const hoNDArray<std::complex<float> > &array, const std::string &name, const int line);
	void print_vector(const std::vector<float> &v);
}

#endif // DEBUG_FUNCTIONS_H
