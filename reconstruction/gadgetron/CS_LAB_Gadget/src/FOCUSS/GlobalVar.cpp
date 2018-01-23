/*	
file name	: 	GlobalVar_FOCUSS.cpp

author		: 	Martin Schwartz	(martin.schwartz@med.uni-tuebingen.de)

version		: 	1.0

date		: 	03.01.2015

description	: 	singleton implementation requires instantiation in cpp file - avoiding linking error
*/

#include "GlobalVar.h"

namespace Gadgetron{
	GlobalVar* GlobalVar::_instance = 0;
}