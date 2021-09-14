 using System;
using Chemistry;
using FlashLFQ;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using IO.ThermoRawFileReader;
using Proteomics.Fragmentation;
using IO.Mgf;
using IO.MzML; 
using MassSpectrometry;
using IO.MzML.Generated;
using MzLibUtil; 

namespace DIAMS1Quant
{
	public static class DIAImport
	{
		public static List<MsDataScan> ImportMGF(string path)
		{
			List<MsDataScan> data = Mgf.LoadAllStaticData(path).GetAllScansList();
			return data; 
		}
		public static List<MsDataScan> ImportMZXML(string path)
		{
			List<MsDataScan> data = Mzml.LoadAllStaticData(path).GetAllScansList();
			
			return data; 
		}
		public static List<MsDataScan> SelectMS1s(List<MsDataScan> allScans)
		{
			List<MsDataScan> ms1Scans = allScans.Where(i => i.MsnOrder == 1).Select(i => i).ToList().OrderBy(i => i.RetentionTime).ToList();
			return ms1Scans; 
			
		}
		public static void FixMissingMS2Fields(List<MsDataScan> allScans, string hcdEnergy, double width)
		{
			DIAScanModifications.UpdateMS2WithMZIsolationWidth(allScans, 1.2);
			DIAScanModifications.UpdateHcdEnergy(allScans, "30");
			DIAScanModifications.RecalculateIsolationRange(allScans);
		}
	}
	public static class DIAMethods
	{
		public static void FillMissingPrecursorScans(List<MsDataScan> listScans)
		{
			int scanIndex = 1;
			for (int i = 0; i < listScans.Count; i++)
			{
				if (listScans[i].MsnOrder != 1)
				{
					listScans[i].SetOneBasedPrecursorScanNumber(scanIndex);
				}
				else
				{
					scanIndex = i + 1; 
				}
			}
		}
	}
	public static class DIAScanModifications
	{
		public static void UpdateMS2WithMZIsolationWidth(List<MsDataScan> allScans, double mzWindow)
		{
			foreach (MsDataScan scan in allScans)
			{
				scan.SetIsolationWidth(mzWindow);
			}
		}
		public static void UpdateHcdEnergy(List<MsDataScan> allScans, string energy)
		{
			foreach(MsDataScan scan in allScans)
			{
				scan.SetHcdEnergy(energy); 
			}
		}
		public static void UpdateIsolationWidth(List<MsDataScan> allScans, double value)
		{
			foreach(MsDataScan scan in allScans)
			{
				scan.SetIsolationWidth(value); 
			}
		}
		public static void RecalculateIsolationRange(List<MsDataScan> allScans)
		{
			foreach(MsDataScan scan in allScans)
			{
				scan.UpdateIsolationRange(); 
			}
		}
		public static void UpdateOneBasedPrecursorScan(List<MsDataScan> combinedScans)
		{
			int currentOneBasedScanNumber = 0;
			foreach (MsDataScan scan in combinedScans)
			{
				if (scan.MsnOrder == 1)
				{
					currentOneBasedScanNumber = scan.OneBasedScanNumber;
				}
				else
				{
					scan.SetOneBasedPrecursorScanNumber(currentOneBasedScanNumber);
				}
			}
		}
		public static void UpdatePrecursorIntensity(List<MsDataScan> combinedScans)
		{
			foreach(MsDataScan scan in combinedScans)
			{
				if(scan.OneBasedPrecursorScanNumber == null)
				{
					break;
				}
				else
				{
					int precursorScanNumber = scan.OneBasedPrecursorScanNumber.Value;
				}
			}
		}

		/* public static List<MsDataScan> UpdateMissingValuesInMS2Spectra(List<MsDataScan> allScans)
		{
			List<MsDataScan> updatedScans = new List<MsDataScan>(); 
			
			foreach(MsDataScan scan in allScans)
			{
				if(scan.MsnOrder != 1)
				{
					MzRange range; 
					if(scan.IsolationRange == null)
					{
						range = MzRange(scan.IsolationMz)
					}

					MsDataScan newScan = new MsDataScan(scan.MassSpectrum, 
						oneBasedPrecursorScanNumber: scan.OneBasedScanNumber, msnOrder: scan.OneBasedScanNumber, 
						isCentroid: scan.IsCentroid, polarity: scan.Polarity, retentionTime: scan.RetentionTime)
				}
				else
				{
					updatedScans.Add(scan); 
				}
			}
		}
		*/
	}
}
