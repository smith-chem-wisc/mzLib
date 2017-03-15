#include "stdafx.h"
#include "Win32Project8.h"
#include "ClassLibrary1.h"
#include < vcclr.h >  

using namespace System::Runtime::InteropServices;

array<ClassLibrary1::PrecursorInfo>^ ClassLibrary1::Class1::runTheMethod(String^ path)
{
	Console::WriteLine(L"Hello World");

	MSFileReaderLib::IXRawfile5Ptr m_Raw = fnWin32Project8();

	pin_ptr<const wchar_t> wch = PtrToStringChars(path);

	m_Raw->Open(wch);


	//TCHAR pth[MAX_PATH];
	//MultiByteToWideChar(CP_ACP, 0, szFileName, -1, (LPWSTR)pth, MAX_PATH);

	//m_Raw->Open((LPWSTR)pth);

	//Set the current controller
	m_Raw->SetCurrentController(0, 1);


	long firstScanNumber = 0;
	m_Raw->GetFirstSpectrumNumber(&firstScanNumber);
	long lastScanNumber = 0;
	m_Raw->GetLastSpectrumNumber(&lastScanNumber);

	Console::WriteLine(L"first scan number : ");
	Console::WriteLine(firstScanNumber);
	Console::WriteLine(L"last scan number : ");
	Console::WriteLine(lastScanNumber);






	//now transform into a .NET object 
	array<PrecursorInfo>^ infos = gcnew array<PrecursorInfo>(lastScanNumber - firstScanNumber+1);


	for (int j = firstScanNumber; j <= lastScanNumber; j++) {

		_variant_t vPrecursorInfos;
		long nPrecursorInfos;
		m_Raw->GetPrecursorInfoFromScanNum(j, &vPrecursorInfos, &nPrecursorInfos);

		//Console::WriteLine("j:" + j);

		BYTE* pData;
		SafeArrayAccessData(vPrecursorInfos.parray, (void**)&pData);

		if (nPrecursorInfos>0)
		{
			MSFileReaderLib::MS_PrecursorInfo precursorInfo;
			memcpy(&precursorInfo,
				pData,
				sizeof(MSFileReaderLib::MS_PrecursorInfo));
			//Console::WriteLine("LnChargeState:" + precursorInfo.nChargeState);
			infos[j - firstScanNumber] = safe_cast<PrecursorInfo>(Marshal::PtrToStructure(IntPtr(&precursorInfo), PrecursorInfo::typeid));

			//Console::WriteLine("LnChargeState:" + infos[j - firstScanNumber]->nChargeState);
		}
		SafeArrayUnaccessData(vPrecursorInfos.parray);
	}

	m_Raw->Close();

	return infos;

}
