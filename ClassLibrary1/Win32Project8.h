// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the WIN32PROJECT8_EXPORTS
// symbol defined on the command line. This symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// WIN32PROJECT8_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef WIN32PROJECT8_EXPORTS
#define WIN32PROJECT8_API __declspec(dllexport)
#else
#define WIN32PROJECT8_API __declspec(dllimport)
#endif

#import "libid:F0C5F3E3-4F2A-443E-A74D-0AABE3237494" rename_namespace("XRawfile") rename("value", "valueEx")
 
using namespace XRawfile;

WIN32PROJECT8_API IXRawfile5Ptr InitializeRawConnection(void);
