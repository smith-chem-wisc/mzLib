using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace InSourceDecay
{
	public interface ScanGroup
	{
		public string MatrixType { get; }
		// The shared functions go here. Could potentially put the matrix functions here. 
		public int FirstScan { get; }
		public int LastScan { get; }
		public double[] IndexMZ { get; }
		public double[] TICValues { get; set; }
		public int ScansToAverage { get; set;  }
		public double AllScansTICVariance { get; set; }
		public void IsValid() { }
	}

	public class ScanCollection<T> : List<T> where T : ScanGroup
	{
		private List<T> ScansList;
		public double AllScanTICVariance { get; private set; }
		public double ScansToAverage { get; private set; }

		public ScanCollection()
		{
			ScansList = new List<T>(); 
		}
		public ScanCollection(List<MsDataScan> allScans)
		{
			ScanCollection<AveragedScan> aveScans = new ScanCollection<AveragedScan>(); 
		}
		public void CalculateAllScansTICVariance(List<MsDataScan> allScans)
		{
			List<double> ticValues = new List<double>();
			double allScansTicVariance = new double(); 
			foreach(MsDataScan scan in allScans)
			{
				ticValues.Add(scan.TotalIonCurrent); 
			}
			allScansTicVariance = MatrixCalculations.Covariance(ticValues.ToArray(), ticValues.ToArray());
			AllScanTICVariance = allScansTicVariance; 
		}
		public void CreateTypeSpecificDataStructure(string matrixType)
		{

		}
		public bool IsAveragedScan()
		{
			if(this.GetType().ToString() == "AveragedScan")
			{
				return true;
			}
			else
			{
				return false; 
			}
		}
		public bool IsPartialCovariance()
		{
			if (this.GetType().ToString() == "PartialCovariance")
			{
				return true;
			}
			else
			{
				return false; 
			}
		}
	}
}
