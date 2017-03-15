#include "stdafx.h"
#include "Win32Project8.h"
#include <iostream>
#include <io.h>
#include <string>

using namespace std;

WIN32PROJECT8_API MSFileReaderLib::IXRawfile5Ptr fnWin32Project8(void)
{
	MSFileReaderLib::IXRawfile5Ptr m_Raw;
	CoInitialize(NULL);
	HRESULT hr = m_Raw.CreateInstance("MSFileReader.XRawfile.1");
	return m_Raw;
}