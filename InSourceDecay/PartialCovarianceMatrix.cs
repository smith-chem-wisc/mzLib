using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InSourceDecay
{
	public class PartialCovarianceMatrix : ScanGroup
	{
		// calculates partial covariance from an AveragedScan object
		// output is a double[IndexMZ.Length, IndexMZ.Length] matrix. 
		public string MatrixType { get; set; }
		public double[,] PCovMatrix { get; set; }
		public int FirstScan { get; set; }
		public int LastScan { get; set; }
		public double[] IndexMZ { get; set; }
		public double[] TICValues { get; set; }
		public int ScansToAverage { get; set; }
		public double AllScansTICVariance { get; set; }

		public PartialCovarianceMatrix(AveragedScan aveScan)
		{
			FirstScan = aveScan.FirstScan;
			LastScan = aveScan.LastScan;
			IndexMZ = aveScan.IndexMZ;
			TICValues = aveScan.TICValues;
			ScansToAverage = aveScan.ScansToAverage;
			AllScansTICVariance = aveScan.AllScansTICVariance;
			CalculatePartialCovarianceMatrix(aveScan);
		}
		public PartialCovarianceMatrix(double[,] partialCovMatrix, int firstScan, int lastScan)
		{
			PCovMatrix = partialCovMatrix;
			FirstScan = firstScan;
			LastScan = lastScan; 
		}
		public void CalculatePartialCovarianceMatrix(AveragedScan aveScan)
		{
			double[,] matrix = aveScan.DataMatrix;
			double[] randomVariable = aveScan.TICValues;
			double allObsVVariance = aveScan.AllScansTICVariance; 

			PCovMatrix = new double[matrix.GetLength(0), matrix.GetLength(0)];
			double[,] covarianceMatrix = MatrixCalculations.Covariance(matrix);
			double[,] partialCovarianceTermMatrix = MatrixCalculations.CalculatePartialCovarianceTerm(matrix, randomVariable, allObsVVariance);
			double[,] partialCovarianceMatrix = MatrixCalculations.Subtract(covarianceMatrix, partialCovarianceTermMatrix);
			PCovMatrix = partialCovarianceMatrix;
		}
		public void IsValid()
		{
			IsValidFirstScan();
			IsValidLastScan();
			IsValidMZIndex();
			IsValidTICVariance();
		}
		public void IsValidMZIndex()
		{
			if (IndexMZ == null)
			{
				throw new ArgumentException("IndexMZ is null");
			}
			foreach (double i in IndexMZ)
			{
				if (i.GetType() != typeof(double))
				{
					throw new ArgumentException("IndexMZ contains object of invalid class");
				}
			}

			int scansOffset = LastScan - FirstScan;
			if (IndexMZ.Length != scansOffset || ScansToAverage != scansOffset)
			{
				throw new ArgumentException("IndexMZ of incorrect length");
			}
		}
		public void IsValidFirstScan()
		{
			if (FirstScan.GetType() != typeof(int))
			{
				throw new ArgumentException("Invalid FirstScan type");
			}
			if (FirstScan % ScansToAverage != 0)
			{
				throw new ArgumentException("Invalid FirstScan Value: Remainder > 0.");
			}
		}
		public void IsValidLastScan()
		{
			if (LastScan.GetType() != typeof(int))
			{
				throw new ArgumentException("Invalid LastScan type");
			}
			if ((LastScan + 1) % 50 != 0)
			{
				throw new ArgumentException("Invalid LastScan Value: Remainder > 0.");
			}
		}
		public void IsValidTICVariance()
		{
			if (AllScansTICVariance.GetType() != typeof(double))
			{
				throw new ArgumentException("Invalid AllScanTICVariance");
			}
		}
		public void TestForNullValuesInPCovMatrix()
		{
			int rows = PCovMatrix.GetLength(0);
			int cols = PCovMatrix.GetLength(1);

			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
				{
					ThrowArgumentExceptionIfNull(PCovMatrix[i, j]);
				}
			}
		}
		public void ThrowArgumentExceptionIfNull(double value)
		{
			if (this.GetType() == null)
			{
				throw new ArgumentException("Value in matrix is null.");
			}
		}

	}
}


