using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry; 

namespace InSourceDecay
{
	public class AveragedScan : ScanGroup
	{
		
		public string MatrixType { get; set; }
		public double[] TICValues { get; set; } 
		// contains all tic values from the scans passed to the object initialization
		public double[,] DataMatrix { get; set; }
		// Array where each column is the intensity values of an MsDataScan and every row position corresponds to an m/z
		public double AllScansTICVariance { get; set; }
		// value passed from AveragedScanCollection
		public int FirstScan { get; set;  }
		// Value of the first scan in the subset 
		public int LastScan { get; set; }
		// Value of the last scan in the subset
		public double[] IndexMZ { get; set; }
		// All unique m/z values in the subset of data passed to object initialization
		public int ScansToAverage { get; set; }
		public AveragedScan(List<MsDataScan> subsetDataScans, int firstScan, double allScanTicVariance, int scansToAverage)
		{
			FirstScan = firstScan + 1;
			LastScan = FirstScan + scansToAverage - 1;
			ScansToAverage = scansToAverage;
			CreateMZIndex(subsetDataScans); // sets the value of the IndexMZ property
			TICValues = AveragedScanCollection.GetTICValues(subsetDataScans).ToArray();
			AllScansTICVariance = allScanTicVariance;
			DataMatrix = CreateDataMatrix(subsetDataScans);

		}
		public void CreateMZIndex(List<MsDataScan> scans)
		{

			List<double> initialIndex = new List<double>();
			List<double[]> xarrays = new List<double[]>();

			foreach (MsDataScan scan in scans)
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
			double[] indexMZ = initialIndex.Distinct(new DoubleEqualityComparer()).ToList().OrderBy(i => i).ToArray();
			IndexMZ = new double[indexMZ.Length];
			IndexMZ = indexMZ;
		}

		public double[,] CreateDataMatrix(List<MsDataScan> subsetScans)
		{
			double[,] dataMatrix = new double[IndexMZ.Length, ScansToAverage];
			for (int k = 0; k < subsetScans.Count; k += ScansToAverage)
			{
				double[] xarray = subsetScans.ElementAt(k).MassSpectrum.XArray;
				double[] yarray = subsetScans.ElementAt(k).MassSpectrum.YArray;

				int i = 0; // iterating through indexMZVector
				int j = 0; // iterating through xarray

				double difference = 0.005;

				while (i != dataMatrix.GetLength(0) && j != xarray.Length)
				{
					if (Math.Abs(IndexMZ[i] - xarray[j]) >= difference)
					{
						dataMatrix[i, k] = yarray[j];
						i++; j++;
					}
					else
					{
						dataMatrix[i, k] = 0;
						i++;
					}
				}
			}
			return dataMatrix; 
		}
		public void IsValid()
		{
			IsValidFirstScan();
			IsValidLastScan();
			IsValidMZIndex();
			IsValidTICVariance();
			IsValidDataMatrix(); 
		}
		public void IsValidMZIndex()
		{
			if(IndexMZ == null)
			{
				throw new ArgumentException("IndexMZ is null"); 
			}
			foreach(double i in IndexMZ)
			{
				if(i.GetType() != typeof(double))
				{
					throw new ArgumentException("IndexMZ contains object of invalid class"); 
				}
			}

			int scansOffset = LastScan - FirstScan; 
			if(IndexMZ.Length != scansOffset || ScansToAverage != scansOffset)
			{
				throw new ArgumentException("IndexMZ of incorrect length"); 
			}
		}
		public void IsValidFirstScan()
		{
			if(FirstScan.GetType() != typeof(int))
			{
				throw new ArgumentException("Invalid FirstScan type"); 
			}
			if(FirstScan % ScansToAverage != 0)
			{
				throw new ArgumentException("Invalid FirstScan Value: Remainder > 0."); 
			}
		}
		public void IsValidLastScan()
		{
			if(LastScan.GetType() != typeof(int))
			{
				throw new ArgumentException("Invalid LastScan type"); 
			}
			if((LastScan + 1) % 50 != 0)
			{
				throw new ArgumentException("Invalid LastScan Value: Remainder > 0."); 
			}
		}
		public void IsValidDataMatrix()
		{
			int rows = DataMatrix.GetLength(0);
			int cols = DataMatrix.GetLength(1);

			if (rows != IndexMZ.Length)
			{
				throw new ArgumentException("DataMatrix row length does not match IndexMZ Length");
			}
			if (cols != ScansToAverage)
			{
				throw new ArgumentException("DataMatrix column length does not match ScansToAverage"); 
			}
			TestForNullValuesInDataMatrix(); 			
		}
		public void IsValidTICVariance()
		{
			if (AllScansTICVariance.GetType() != typeof(double))
			{
				throw new ArgumentException("Invalid AllScanTICVariance"); 
			}
		}
		public void TestForNullValuesInDataMatrix()
		{
			int rows = DataMatrix.GetLength(0);
			int cols = DataMatrix.GetLength(1); 

			for(int i = 0; i < rows; i++)
			{
				for(int j = 0; j < cols; j++)
				{
					ThrowArgumentExceptionIfNull(DataMatrix[i,j]); 
				}
			}
		}
		public void ThrowArgumentExceptionIfNull(double value)
		{
			if(this.GetType() == null)
			{
				throw new ArgumentException("Value in matrix is null."); 
			}
		}

	}
}
