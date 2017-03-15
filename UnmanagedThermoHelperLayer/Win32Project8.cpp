#include "stdafx.h"
#include "Win32Project8.h"
#include <iostream>
#include <io.h>
#include <string>

using namespace std;

WIN32PROJECT8_API IXRawfile5Ptr InitializeRawConnection(void)
{
	IXRawfile5Ptr m_Raw;
	CoInitialize(NULL);
	HRESULT hr = m_Raw.CreateInstance("MSFileReader.XRawfile.1");
	return m_Raw;
}