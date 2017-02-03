#include "CS_GlobalVar.h"

// singleton implementation requires instantiation in cpp file - avoiding linking error

namespace Gadgetron{
	CS_GlobalVar* CS_GlobalVar::_instance = 0;
}