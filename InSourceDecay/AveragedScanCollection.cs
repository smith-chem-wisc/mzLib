using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry; 

namespace InSourceDecay
{
	public class AveragedScanCollection
	{
		/* This class contains the data and functions to perform a 2D-PC-MS workflow: 
		 * 1) Calculate the all scan TIC variance
		 * 2) Create an AveragedScan object that contains all the information required to calculate
		 * the partial covariance for a set of scans. 
		 * 3) Calculate the partial covariance matrix and return a List<PartialCovarianceMatrix> for further analysis
		 */
		public double TICVariance { get; set; }
		public double[] TICValues { get; set; }
		public List<AveragedScan> CombinedScans { get; set; }
		public int ScansToAverage { get; set; }
		public List<PartialCovarianceMatrix> PartialCovarianceList { get; set; } 

		public AveragedScanCollection(List<MsDataScan> allScans, int scansToAverage)
		{
			TICValues = GetTICValues(allScans).ToArray();
			TICVariance = MatrixCalculations.Covariance(TICValues, TICValues);
			ScansToAverage = scansToAverage;
			CombinedScans = new List<AveragedScan>();
			PartialCovarianceList = new List<PartialCovarianceMatrix>(); 
		}




		public static List<double> GetTICValues(List<MsDataScan> msScans)
		{
			List<double> TICValues = new List<double>();
			foreach (MsDataScan scan in msScans)
			{
				TICValues.Add(scan.TotalIonCurrent);
			}
			double[] TICValuesArray = TICValues.ToArray();
			return TICValuesArray.ToList();
		}
		public static List<double> CombineMZValuesToIndex(List<MsDataScan> allScans)
		{
			List<double> initialIndex = new List<double>();
			List<double[]> xarrays = new List<double[]>();

			foreach (MsDataScan scan in allScans)
			{
				xarrays.Add(scan.MassSpectrum.XArray);
			}
			foreach (double[] vector in xarrays)
			{
				for (int i = 0; i < vector.Length; i++)
				{
					initialIndex.Add(vector[i]);
				}
			}
			return initialIndex.Distinct(new DoubleEqualityComparer()).ToList().OrderBy(i => i).ToList();
		}
		public AveragedScan CreateAveragedScan(List<MsDataScan> scans, int scanStart)
		{
			// this scanStart is the 0-based index value, not the 1-based index used to denote scan position in time. 
			List<MsDataScan> scansSubset = scans.Skip(scanStart).Take(ScansToAverage).ToList();
			AveragedScan averagedScan = new AveragedScan(scansSubset, scanStart, TICVariance, ScansToAverage);
			return averagedScan;
		}
		public void CreateAverageScanList(List<MsDataScan> allScans)
		{
			// will probably need to parallelize this code well in order to get it to work in a timely manner
			
			for (int i = 0; i < allScans.Count; i += ScansToAverage)
			{
				CombinedScans.Add(CreateAveragedScan(allScans, i));
			}
		}
		public void CalculatePartialCovarianceMatrices()
		{
			// foreach matrix, calculate the covariance
			// Calculate matrix covariance
			// Calculate the partial covariance term
			// Covariance - partial covariance term
			foreach (AveragedScan aveScan in CombinedScans)
			{
				PartialCovarianceMatrix pCovMatrix = new PartialCovarianceMatrix(aveScan);
				PartialCovarianceList.Add(pCovMatrix); 
			}
		}
		public void CalculatePartialCovarianceMatricesParallel()
		{
			Parallel.ForEach(CombinedScans, scan =>
			{
				PartialCovarianceMatrix pCovMatrix = new PartialCovarianceMatrix(scan);
				PartialCovarianceList.Add(pCovMatrix);
			}); 
		}
		public class AveragedScanCollectionEnumerator
		{
			int nIndex;
			AveragedScanCollection averageScanCollection; 
			public AveragedScanCollectionEnumerator(AveragedScanCollection coll)
			{
				averageScanCollection = coll;
				nIndex = -1; 
			}
			public bool MoveNext()
			{
				nIndex++;
				return (nIndex < averageScanCollection.CombinedScans.Count); 
			}
			public AveragedScan Current => averageScanCollection.CombinedScans[nIndex]; 
		}
	}
}
