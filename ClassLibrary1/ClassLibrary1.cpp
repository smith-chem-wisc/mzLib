#include "stdafx.h"
#include "Win32Project8.h"
#include "ClassLibrary1.h"
#include < vcclr.h >  

using namespace System::Runtime::InteropServices;

array<ManagedThermoHelperLayer::PrecursorInfo>^ ManagedThermoHelperLayer::HelperClass::GetAllPrecursorInfos(String^ path)
{
	IXRawfile5Ptr m_Raw = InitializeRawConnection();

	pin_ptr<const wchar_t> wch = PtrToStringChars(path);

	m_Raw->Open(wch);

	m_Raw->SetCurrentController(0, 1);

	long firstScanNumber = 0;
	m_Raw->GetFirstSpectrumNumber(&firstScanNumber);
	long lastScanNumber = 0;
	m_Raw->GetLastSpectrumNumber(&lastScanNumber);

	//now transform into a .NET object 
	array<PrecursorInfo>^ infos = gcnew array<PrecursorInfo>(lastScanNumber - firstScanNumber+1);

	for (int j = firstScanNumber; j <= lastScanNumber; j++) {

		_variant_t vPrecursorInfos;
		long nPrecursorInfos;
		m_Raw->GetPrecursorInfoFromScanNum(j, &vPrecursorInfos, &nPrecursorInfos);
		
		BYTE* pData;
		SafeArrayAccessData(vPrecursorInfos.parray, (void**)&pData);

		if (nPrecursorInfos>0)
		{
			MS_PrecursorInfo precursorInfo;

			// Always grab the FIRST precursor info to be the actual precursor info
			memcpy(&precursorInfo,
				pData,
				sizeof(MS_PrecursorInfo));
			infos[j - firstScanNumber] = safe_cast<PrecursorInfo>(Marshal::PtrToStructure(IntPtr(&precursorInfo), PrecursorInfo::typeid));
		}
		SafeArrayUnaccessData(vPrecursorInfos.parray);
	}

	m_Raw->Close();

	return infos;

}
