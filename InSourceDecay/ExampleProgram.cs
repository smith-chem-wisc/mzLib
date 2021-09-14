using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System.IO;
using IO.ThermoRawFileReader;
using Proteomics.Fragmentation;
using IO.Mgf;
using InSourceDecay;

namespace InSourceDecay
{
	class ExampleProgram
	{
		static void Main(string[] args)
		{
			// Data Import
			string DataFilePath = Path.GetFullPath(@"Smarter_DDA_Top10_SID75_RES_30k.mgf");
			List<MsDataScan> allScans = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();
			// Take a subset of scans to reduce processing time: 
			List<MsDataScan> subsetScans = allScans.Skip(1000).Take(50).ToList();

			// Create an AveragedScanCollection object
			AveragedScanCollection scanCollection = new AveragedScanCollection(subsetScans, 5);
			// Create the list of averaged scans within scanCollection
			scanCollection.CreateAverageScanList(subsetScans); // could probably use a speed-up
			// Create the partial covariance matrices
			scanCollection.CalculatePartialCovarianceMatricesParallel(); 
			// Need a way to easily interrogate the results of the partial covariance matrix. Not sure how that's 
			// going to work yet. Potentially a tuple? 
		}
	}
}
