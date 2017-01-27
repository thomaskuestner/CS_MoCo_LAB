// singleton class implementation for global variables
// ref to http://de.wikibooks.org/wiki/C%2B%2B-Programmierung:_Entwurfsmuster:_Singleton

#ifndef CS_GLOBALVAR_H
#define CS_GLOBALVAR_H

#include "Gadget.h"
#include "NDArray.h"
#include "GadgetMRIHeaders.h"
#include "ismrmrd.h"

namespace Gadgetron{

	class CS_GlobalVar{
		public:
			static CS_GlobalVar* instance(){
				static CGuard g;
				if (!_instance)
					_instance = new CS_GlobalVar();
				return _instance;
			}

		// acquisition header
		std::vector<GadgetContainerMessage<ISMRMRD::AcquisitionHeader>*> AcqVec_;

		// image header
		std::vector<GadgetContainerMessage<ISMRMRD::ImageHeader>*> ImgHeadVec_;

		private:
			static CS_GlobalVar* _instance;
			CS_GlobalVar() {}
			CS_GlobalVar(const CS_GlobalVar&);
			~CS_GlobalVar(){}

			class CGuard{
				public:
					~CGuard()
					{
						if( NULL != CS_GlobalVar::_instance)
						{
							delete CS_GlobalVar::_instance;
							CS_GlobalVar::_instance = NULL;
						}
					}
			};
			friend class CGuard;
	};	
}

#endif  //CS_GLOBALVAR_H