using System;
using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using IO.ThermoRawFileReader;
using Proteomics.Fragmentation;
using MathNet;
using System.Threading.Tasks;

namespace InSourceDecay
{
	// Used for devleopment of the workflow. Not actually used in the program. 
	public class PartialCovarianceCalculations
	{
		// Data pre-processing
		// interating with intial objects
		// Creates the data matrix with List<double> and List<MsDataScan>
		public double[,] CreateDataMatrix(List<double> indexMZvector, List<MsDataScan> allScans)
		{
			double[,] dataMatrix = new double[indexMZvector.Count, allScans.Count];

			for (int k = 0; k < allScans.Count; k++)
			{

				double[] xarray = allScans.ElementAt(k).MassSpectrum.XArray;
				double[] yarray = allScans.ElementAt(k).MassSpectrum.YArray;

				int i = 0; // iterating through indexMZVector
				int j = 0; // iterating through xarray

				double difference = 0.005;

				while (i != dataMatrix.GetLength(0) && j != xarray.Length)
				{
					if (Math.Abs(indexMZvector[i] - xarray[i]) >= difference)
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
		// Creates a data matrix, using only List<MsDataScan> as input
		public double[,] CreateDataMatrix(List<MsDataScan> allScans)
		{
			PartialCovarianceCalculations partialCov = new PartialCovarianceCalculations();
			List<double> indexMZvector = partialCov.CombineMZValuesToIndex(allScans);

			double[,] dataMatrix = new double[indexMZvector.Count, allScans.Count];

			for (int k = 0; k < allScans.Count; k++)
			{

				double[] xarray = allScans.ElementAt(k).MassSpectrum.XArray;
				double[] yarray = allScans.ElementAt(k).MassSpectrum.YArray;

				int i = 0; // iterating through indexMZVector
				int j = 0; // iterating through xarray

				double difference = 0.005;

				while (i != dataMatrix.GetLength(0) && j != xarray.Length)
				{
					if (Math.Abs(indexMZvector[i] - xarray[j]) >= difference)
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
		// fetches all the x arrays froma List<MsDataScan> and puts it into a list<double[]>
		public static List<double[]> AccessXArrays(List<MsDataScan> scanslist, List<double[]> xarrays)
		{

			foreach (MsDataScan scan in scanslist)
			{
				xarrays.Add(scan.MassSpectrum.XArray);
			}
			return xarrays;
		}
		// fetches all the xarrays and puts it into a List<double[]>. Doesn't require the output list as a variable. 
		public static List<double[]> AccessXarrays(List<MsDataScan> scansList)
		{
			List<double[]> xarrays = new List<double[]>(); 
			foreach(MsDataScan scan in scansList)
			{
				xarrays.Add(scan.MassSpectrum.XArray); 
			}
			return xarrays; 
		}
		// fetches all the Yarrays from a List<MsDataScan> 
		public static List<double[]> FetchIntensityArrays(List<MsDataScan> allScans)
		{
			List<double[]> intensityArrays = new List<double[]>();
			foreach (MsDataScan scan in allScans)
			{
				intensityArrays.Add(scan.MassSpectrum.YArray);
			}
			return intensityArrays;
		}
		public List<double> GetTICValues(List<MsDataScan> msScans)
		{
			List<double> TICValues = new List<double>();
			foreach (MsDataScan scan in msScans)
			{
				TICValues.Add(scan.TotalIonCurrent);
			}
			double[] TICValuesArray = TICValues.ToArray();
			return TICValuesArray.ToList();
		}
		// Creates a single List<double> of all m/z values from all scans included in List<double[]> format
		public static List<double> CombineMZValuesToIndex(List<double[]> combinedXArrays)
		{
			List<double> initialIndex = new List<double>();
			foreach (double[] xarray in combinedXArrays)
			{
				for (int i = 0; i < xarray.Length; i++)
				{
					initialIndex.Add(xarray[i]);
				}

			}
			return initialIndex.Distinct(new DoubleEqualityComparer()).ToList().OrderBy(i => i).ToList();
		}
		// Creates a single List<double> of all m/z values from all scans in List<MsDataScan> format. 
		public List<double> CombineMZValuesToIndex(List<MsDataScan> allScans)
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
		public static List<double[]> ReturnValuesFromListOfArray(List<double[]> ListOfVectors, List<double> index)
		{
			List<double[]> outputList = new List<double[]>();
			for (int i = 0; i < index.Count; i++)
			{
				double[] tempCollectionArray = ListOfVectors.Select(a => a[i]).ToArray();
				outputList.Add(tempCollectionArray);
			}
			return outputList;
		}

		// Creates a list of 2D arrays where the columns are each scan and the column position is m/z value
		public List<double[,]> AverageScans(double[,] dataMatrix, int scansToAverage)
		{
			int numberOfScans = dataMatrix.GetLength(1);
			int lastIterationValue = numberOfScans - (scansToAverage - 1);
			List<double[,]> averagedScansCollection = new List<double[,]>();
			for (int i = 0; i < dataMatrix.GetLength(0); i++)
			{
				double[,] tempAveragedScan = new double[dataMatrix.GetLength(0), numberOfScans];
				for (int j = 0; j < scansToAverage; j++)
				{
					tempAveragedScan[i, j] = dataMatrix[i, j];
				}
				averagedScansCollection.Add(tempAveragedScan);
			}

			return averagedScansCollection;
		}

		public static double[] UpdateSpectraWithZeroes(List<double> indexMZvector, MsDataScan scan)
		{
			double[] allValuesVector = new double[indexMZvector.Count];
			double[] xarray = scan.MassSpectrum.XArray;
			double[] yarray = scan.MassSpectrum.YArray;

			int i = 0; // iterating through indexMZVector
			int j = 0; // iterating through xarray

			double difference = 0.005;
			while (i != allValuesVector.Length && j != xarray.Length)
			{
				if (Math.Abs(allValuesVector[i] - xarray[i]) >= difference)
				{
					allValuesVector[i] = yarray[j];
					i++;
					j++;
				}
				else
				{
					allValuesVector[i] = 0;
					i++;
				}
			}
			return allValuesVector;
		}

		// Covariance Calculations


		

		// Correlation
		public void CalcCorrelationScore()
		{

		}
		public double Correlation(double[] array1, double[] array2)
		{
			double avg1 = array1.Average();
			double avg2 = array2.Average();

			double sum1 = array1.Zip(array2, (x1, y1) => (x1 - avg1) * (y1 - avg2)).Sum();

			double sumSqr1 = array1.Sum(x => Math.Pow((x - avg1), 2.0));
			double sumSqr2 = array2.Sum(y => Math.Pow((y - avg2), 2.0));

			double result = sum1 / Math.Sqrt(sumSqr1 * sumSqr2);
			return result;
		}
		public static void Write2DArrayToTestContext(double[,] array2D)
		{
			for (int i = 0; i < array2D.GetLength(0); i++)
			{
				for (int j = 0; j < array2D.GetLength(0); j++)
				{
					TestContext.WriteLine(array2D[i, j].ToString());
				}
			}
		}

	}

	public class PartialCovarianceHelpers
	{
		public static double GetElementAtPosition(List<double[,]> list, int element, int dim, int index)
		{
			double result = (double)list.ElementAt(element).GetValue(dim, index);
			return (result);
		}
	}

	public class MultiScanCovariance
	{
		public List<double> CombinedMassToCharge { get; set; }
		public List<double> Intensities1 { get; set; }
		public List<double> Intensities2 { get; set; }

		MultiScanCovariance(List<MsDataScan> msScans)
		{
		}
	}
	class DoubleEqualityComparer : IEqualityComparer<double>
	{
		public bool Equals(double d1, double d2)
		{
			double difference = 0.005; 
			if(Math.Abs(d1 - d2) <= difference)
			{
				return true; 
			}else if (Math.Abs(d1-d2) > difference)
			{
				return false;
			}
			else
			{
				return false; 
			}
		}
		public int GetHashCode(double d1)
		{
			return (int) Math.Round(d1); 
		}
	}
}
