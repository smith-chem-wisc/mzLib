using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MassSpectrometry; 
using DIAMS1Quant;
using NUnit.Framework;
using IO.Mgf;
using IO.MzML;
using IO.MzML.Generated; 
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using System.Globalization;
using System.Security.Cryptography;
using System.Text.RegularExpressions;
using System.Xml.Serialization;
using System.ComponentModel;


namespace Test
{
	public class DIAImportTests
	{
		public void PrintMsDataScanProperties(MsDataScan scan)
		{
			foreach (PropertyDescriptor descriptor in TypeDescriptor.GetProperties(scan))
			{
				string name = descriptor.Name;
				object value = descriptor.GetValue(scan);
				Console.WriteLine("{0}; {1}", name, value);
			}
		}
		[Test]
		public void TestDataProcessing()
		{
			string filePathmzML = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
				"20210211_FT_nDIA_1.mzML");
			string filePathMgf = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles",
				"MS2_data.mgf");

			List<MsDataScan> ms1Scans = DIAImport.SelectMS1s(DIAImport.ImportMZXML(filePathmzML));
			List<MsDataScan> ms2Scans = DIAImport.ImportMGF(filePathMgf);
			DIAImport.FixMissingMS2Fields(ms2Scans, "30", 1.2); 

			List<MsDataScan> combinedScans = ms1Scans.Concat(ms2Scans).ToList().OrderBy(i => i.RetentionTime).ToList();

			for (int i = 0; i < combinedScans.Count; i++)
			{
				combinedScans[i].SetOneBasedScanNumber(i + 1);
			}

			int currentOneBasedScanNumber = 0; 
			foreach(MsDataScan scan in combinedScans)
			{
				if(scan.MsnOrder == 1)
				{
					currentOneBasedScanNumber = scan.OneBasedScanNumber;
				}
				else
				{
					scan.SetOneBasedPrecursorScanNumber(currentOneBasedScanNumber); 
				}
			}

			foreach(MsDataScan scan in combinedScans)
			{
				TestContext.Write(scan.MsnOrder.ToString() + "; " + scan.RetentionTime.ToString() + "\n"); 

			}

			/*
			List<MsDataScan> subsetScans = combinedScans.Skip(5500).Take(3000).ToList();
			foreach (MsDataScan scan in subsetScans)
			{
				PrintMsDataScanProperties(scan);
			}


			// put scans into an array instead of a list 
			MsDataScan[] combinedScansArray = combinedScans.ToArray();
			/*

						// get my file reader
			var myFileReader = new MyFileManager(true);
			// read in one of the original file paths 
			var commonParameters = new CommonParameters();
			MsDataFile outputFile = myFileReader.LoadFile(filePathmzML, commonParameters);
			MsDataFile combinedScansOutputFile = new MsDataFile(combinedScansArray, outputFile.SourceFile);
			
			FakeMsDataFile fakeFile = new FakeMsDataFile(subsetScans.ToArray()); 

			// write file
			// string outputFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "truncated_outputMZML.mzML");
			// MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(fakeFile, outputFilePath, false); 
			*/  
		}
	}
}