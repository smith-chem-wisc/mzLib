// ClassLibrary1.h

#pragma once

using namespace System;

namespace ManagedThermoHelperLayer {

	public value struct PrecursorInfo
	{
		double dIsolationMass;
		double dMonoIsoMass;
		long nChargeState;
		long nScanNumber;
	};

	public ref class HelperClass
	{
		// TODO: Add your methods for this class here.
	public:
		array<PrecursorInfo>^ GetAllPrecursorInfos(String^ path);
	};
}
