//CS_LAB_export.h

#ifndef CS_LAB_EXPORT_H_
#define CS_LAB_EXPORT_H_

#if defined (WIN32)
#if defined (CS_LAB_EXPORTS)
#define EXPORTCSLAB __declspec(dllexport)
#else
#define EXPORTCSLAB __declspec(dllimport)
#endif
#else
#define EXPORTCSLAB
#endif

#endif /* CS_LAB_EXPORT_H_ */