// ClassLibrary1.h

#pragma once

using namespace System;

namespace ClassLibrary1 {

	public value struct PrecursorInfo
	{
		double dIsolationMass;
		double dMonoIsoMass;
		long nChargeState;
		long nScanNumber;
	};

	public ref class Class1
	{
		// TODO: Add your methods for this class here.
	public:
		array<PrecursorInfo>^ runTheMethod(String^ path);
	};
}
