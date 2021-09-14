using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using IO.ThermoRawFileReader;
using Proteomics.Fragmentation;
using IO.Mgf;
using InSourceDecay; 

namespace Test
{
	[TestFixture]
	class TestPartialCovarianceCalculations
	{
		public static void PrintMatrix(double[,] matrix, bool testContext)
		{
			for (int i = 0; i < matrix.GetLength(0); i++)
			{
				for (int j = 0; j < matrix.GetLength(1); j++)
				{
					if (testContext != true)
					{
						Console.Write(matrix[i, j].ToString() + " ");
					}
					else
					{
						TestContext.Write(matrix[i, j].ToString() + " ");
					}
				}
				if (testContext != true)
				{
					Console.Write("\n");
				}
				else
				{
					TestContext.Write("\n");
				}
			}
		}
		public static void CheckMatrixEqualityInUnitTest(double[,] result,
			double[,] expected,
			double delta = 0.0001)
		{
			for (int i = 0; i < result.GetLength(0); i++)
			{
				for (int j = 0; j < result.GetLength(1); j++)
				{
					Assert.AreEqual(expected[i, j], result[i, j], delta);
				}
			}
		}
		[Test]
		public void top200peak()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "ISD_test.raw");
			var filterParams = new FilteringParams(200, 0.01, 0, 1, false, true, true);
			var data = ThermoRawFileReader.LoadAllStaticData(DataFilePath, filterParams);
			var testScan = data.GetOneBasedScan(3769);

			double[] peaks = testScan.MassSpectrum.XArray;

			string[] mzMatrix = new string[peaks.Length];
			for (int i = 0; i < peaks.Length; i++)
			{
				double[] temp = new double[peaks.Length];
				for (int j = 0; j < peaks.Length; j++)
				{
					temp[j] = peaks[i] - peaks[j];
				}
				mzMatrix[i] = String.Join(",", temp);
			}

			File.WriteAllLines("C:/Users/Austin/Desktop/output.csv", mzMatrix);

		}
		[Test]
		public void TheoreticalSpectra()
		{
			// TO DO: 
			// Take a protein sequence and digest it
			// Add modifications corresponding to N- or C-terminal mass shifts
			// Combine to make a FASTA
			string proteinSequenceThioRedoxin = "TTFNIQDGPDFQDRVVNSETPVVVDFHAQWCGPCKILGPRLEKMVAKQHGKVVMAKVDIDDHTDLAIEYEVSAVPTVLAMKNGDVVDKFVGIKDEDQLEAFLKKLIG";
			Protein Thioredoxin = new Protein(proteinSequenceThioRedoxin, null);
			DigestionParams digestionParams = new DigestionParams(protease: "non-specific",
				maxMissedCleavages: 35,
				minPeptideLength: 7,
				maxPeptideLength: 34);

			List<Modification> variableModifications = new List<Modification>();
			var digestedThio = Thioredoxin.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();
			List<Product> products = new List<Product>();
			foreach (PeptideWithSetModifications pep in digestedThio)
			{
				pep.Fragment(DissociationType.EThcD, FragmentationTerminus.Both, products);
			}
		}
		[Test]
		public void TestCalcMultiScanCovariance()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "ISD_test.raw");
			// var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			var data = ThermoRawFileReader.LoadAllStaticData(DataFilePath);
			List<MsDataScan> allScans = data.GetAllScansList();

			allScans.ElementAt(1).MassSpectrum.CopyTo2DArray();
			double scanTIC = allScans.ElementAt(1).TotalIonCurrent;

			for (int i = 0; i < 100; i++)
			{
				TestContext.WriteLine(allScans.ElementAt(i).MassSpectrum.XArray.ToString());
			}
		}
		/*
		[Test]
		public void TestCreateMzIndex()
		{
			// need: 
			// List<MsDataScan> scansList
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false); 
			var data = Mgf.LoadAllStaticData(DataFilePath, filterParams).GetAllScansList();

			List<double> indexVector =  PartialCovarianceCalculations.CreateMZIndexVector(data);
			indexVector.ForEach(i => Console.WriteLine(string.Join("\n", i))); 
		}
		*/
		/*[Test]
		public void TestCreateMzindexTolerance()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			var data = Mgf.LoadAllStaticData(DataFilePath, filterParams).GetAllScansList();
			
			List<double> indexVector = PartialCovarianceCalculations.CreateMZIndexWithTolerance(data);
			indexVector.ForEach(i => Console.WriteLine(string.Join("\n", i)));
		}*/
		[Test]
		public void TestCombineMZValuesToIndex()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			var data = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();
			List<double[]> xarrays = new List<double[]>();
			List<double[]> listXarrays = PartialCovarianceCalculations.AccessXArrays(data, xarrays);

			List<double> output = PartialCovarianceCalculations.CombineMZValuesToIndex(listXarrays);
			foreach(double i in output)
			{
				TestContext.WriteLine(string.Join(", ", i, "\n")); 
			} 
		}
		[Test]
		public void TestCovarianceCalc()
		{
			double[] array1 = { 0.3, 0.2, 0.1, -0.4, 0.8 };
			double[] array2 = { 0.5, 0.8, -0.7, -0.25, 0.6 };

			double result = MatrixCalculations.Covariance(array1, array2);
			double expected = 0.1575;
			Assert.AreEqual(expected, result, 0.0001); 
		}

		[Test]
		public void TestCorrelation()
		{
			PartialCovarianceCalculations partialCov = new PartialCovarianceCalculations();
			double[] array1 = { 0.3, 0.2, 0.1, -0.4, 0.8 };
			double[] array2 = { 0.5, 0.8, -0.7, -0.25, 0.6 };

			double result = partialCov.Correlation(array1, array2);
			double expected = 0.5750416;
			Assert.AreEqual(expected, result, 0.00001);
		}
		[Test]
		public void PositiveDefinite()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			var data = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();
			List<double[]> xarrays = new List<double[]>();
			List<double[]> listXarrays = PartialCovarianceCalculations.AccessXArrays(data, xarrays);

			List<double> indexXarray = PartialCovarianceCalculations.CombineMZValuesToIndex(listXarrays);

			double[] outputYarray = PartialCovarianceCalculations.UpdateSpectraWithZeroes(indexXarray, data.ElementAt(3));
			List<double[]> listYarrays = new List<double[]>();
			foreach (MsDataScan scan in data)
			{
				listYarrays.Add(PartialCovarianceCalculations.UpdateSpectraWithZeroes(indexXarray, scan));
			}
			List<double[]> transposedYarrays = PartialCovarianceCalculations.ReturnValuesFromListOfArray(listYarrays.Take(10).ToList(), indexXarray);
			transposedYarrays.ToArray(); 

		}
		[Test]
		public void TestCreateDataMatrix()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			var data = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();

			PartialCovarianceCalculations partialCov = new PartialCovarianceCalculations();
			double[,] dataMatrix = partialCov.CreateDataMatrix(data.Skip(1500).Take(1000).ToList());
		}

		[Test]
		public void TestScanAveragingCorrelationMatrices()
		{
			// I guess I should break up each set of scans into an averaged scan. Probably should be dependent on the "first protein elution or something" 

		}
		[Test]
		public void TestCovarianceWithMatrixInput()
		{
			double[,] testDoubleMatrix = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 }, { 10, 11, 12 } };
			PartialCovarianceCalculations partialCov = new PartialCovarianceCalculations();
			double[,] result = MatrixCalculations.Covariance(testDoubleMatrix);
			double[,] expected = { { 15, 15, 15},
			{ 15,15,15},
			{ 15,15,15} }; 

			for(int i = 0; i < result.GetLength(0); i++)
			{
				for(int j = 0; j < result.GetLength(1); j++)
				{
					Assert.AreEqual(expected[i, j], result[i, j], 0.001); 
				}
			}
		}
		[Test]
		public void TestTranspose()
		{
			double[,] testMatrix = { { 1, 0, 1 }, { 0, 1, 0 } };
			double[,] expected = { { 1, 0 }, { 0, 1 }, { 1, 0 } };
			double[,] result = MatrixCalculations.Transpose(testMatrix);

			Assert.AreEqual(expected[1, 1], result[1, 1]); 
		}
		[Test]
		public void TestFillMatrixOfOnes()
		{
			double[,] result = MatrixCalculations.FillMatrixOf1s(2, 2);
			Assert.AreEqual(1, result[1, 1]); 
		}
		[Test]
		public void TestMatrixMultiply()
		{
			double[,] testArray1 = { { 1, 2, 3 }, { 4, 5, 6 } };
			double[,] testArray2 = { { 1, 2 }, { 3, 4 }, { 5, 6 } };

			double[,] result = MatrixCalculations.Multiply(testArray1, testArray2);
			double[,] expected = { { 22, 28 }, { 49, 64 } };

			for (int i = 0; i < result.GetLength(0); i++)
			{
				for (int j = 0; j < result.GetLength(1); j++)
				{
					Assert.AreEqual(expected[i, j], result[i, j], 0.001); 
				}
			}
		}
		[Test]
		public void TestCovarianceMatrixXVector()
		{
			double[,] testArray1 = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
			double[] testVector = { 1.5, 2.5, 3.5 };
			
			double[] result = MatrixCalculations.Covariance(testArray1, testVector);
			double[] expected = { 1, 1, 1 }; 
			for (int i = 0; i < result.GetLength(0); i++)
			{
				Assert.AreEqual(expected[i], result[i]); 
			}
		}
		[Test]
		public void TestTwoVectorCovarianceToMatrix()
		{
			double[] testArray1 = { 1, 2, 3, 6, -1, -5 };
			double[] testArray2 = { -3, 4, -2, -1, 0, 12 };
			double[,] expected = {
				{ -3, -6, -9, -18, 3, 15},
				{ 4, 8, 12, 24, -4, -20 },
				{ -2, -4, -6, -12, 2, 10},
				{ -1, -2, -3, -6, 1, 5 },
				{ 0, 0, 0, 0, 0 , 0 },
				{12, 24, 36, 72, -12, -60 }
			}; // messed up row/column orientation when creating 2D arrays in C#. 
			// Matrix needs to be transposed to be correct. 
			double[,] result = MatrixCalculations.Multiply(testArray1, testArray2);
			CheckMatrixEqualityInUnitTest(result, 
				MatrixCalculations.Transpose(expected)); 
		}
		[Test]
		public void TestAveragedScanCollectionInitialization()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			List<MsDataScan> allScans = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();
			/* TestContext.Write(scanCollection.CombinedScans.Count().ToString());
			TestContext.Write("\n"); 
			TestContext.Write(scanCollection.CombinedScans.ElementAt(1).FirstScan.ToString());
			TestContext.Write("\n");
			TestContext.Write(scanCollection.CombinedScans.ElementAt(1).LastScan.ToString());
			TestContext.Write("\n");
			TestContext.Write(string.Join(", ", scanCollection.CombinedScans.ElementAt(1).TICValues));
			
			foreach(AveragedScan aveScan in scanCollection.CombinedScans)
			{
				TestContext.Write(aveScan.FirstScan.ToString());
				TestContext.Write("; "); 
				TestContext.Write(aveScan.LastScan.ToString());
				TestContext.Write("\n"); 
			}
			*/
			// scanCollection.CalculatePartialCovarianceMatrices();
			AveragedScanCollection scanCollection = new AveragedScanCollection(allScans, 50); 
		}

		[Test]
		public void TestCreateAverageScanList()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			List<MsDataScan> allScans = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();
			List<MsDataScan> subsetScans = allScans.Skip(1000).Take(50).ToList(); 
			AveragedScanCollection scanCollection = new AveragedScanCollection(subsetScans, 5);
			scanCollection.CreateAverageScanList(subsetScans);
			TestContext.WriteLine(scanCollection.CombinedScans.Count.ToString()); 
		}

		[Test]
		public void TestPartialCovarianceMatrix()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			List<MsDataScan> allScans = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();
			List<MsDataScan> testSubsetScans = allScans.Skip(1500).Take(100).ToList(); 
			AveragedScanCollection scanCollection = new AveragedScanCollection(testSubsetScans, 5);
			scanCollection.CreateAverageScanList(testSubsetScans);
			var testScan = scanCollection.CombinedScans.ElementAt(3);
			PartialCovarianceMatrix result = new PartialCovarianceMatrix(testScan); 
		}

		[Test]
		public void TestPartialCovarianceListCreation()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			List<MsDataScan> allScans = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList().Skip(1000).Take(100).ToList();
			AveragedScanCollection scanCollection = new AveragedScanCollection(allScans, 10);
			scanCollection.CreateAverageScanList(allScans); 
			scanCollection.CalculatePartialCovarianceMatricesParallel();
			TestContext.Write(scanCollection.PartialCovarianceList[8].PCovMatrix.GetLength(0).ToString() +
				scanCollection.PartialCovarianceList[8].PCovMatrix.GetLength(1).ToString()); 
		}

		[Test]
		public void TestCreateAverageScan()
		{
			string DataFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Smarter_DDA_Top10_SID75_RES_30k.mgf");
			var filterParams = new FilteringParams(200, 0.0001, 0, 1, false, false, false);
			List<MsDataScan> allScans = Mgf.LoadAllStaticData(DataFilePath).GetAllScansList();
			List<MsDataScan> subsetScans = allScans.Skip(1000).Take(5).ToList(); 
			AveragedScan testAveScan = new AveragedScan(subsetScans, 0, 6e6, 5);

		}

		[Test]
		public void TestMultiplyTwoVectors()
		{
			double[] vector1 = { 1, 2, 3, 4, 5 };
			double[] vector2 = { 8, 7, 6, 5, 4 };
			double[,] result = MatrixCalculations.Multiply(vector1, vector2);
			PrintMatrix(result, true); 
		}

		[Test]
		public void TestCalculatePartialCovarianceTerm()
		{
			double[,] testMatrix = { { 8, 7, 6, 5, 4 }, { 16, 14, 12, 10, 8 }, { 24, 21, 18, 15, 12 }, { 32, 28, 24, 20, 16 }, { 40, 35, 30, 25, 20 } };
			double[] testRandomVariable = { 4, 10, 12, 16, 20 };
			double varianceRandomVariable = MatrixCalculations.Covariance(testRandomVariable, testRandomVariable);
			double[,] result = MatrixCalculations.CalculatePartialCovarianceTerm(testMatrix, testRandomVariable, varianceRandomVariable);
			PrintMatrix(result, true); 
		}

		[Test] 
		public void TestCalcualtePartialCovarianceMatrix()
		{
			double[,] testMatrix = { { 8, 7, 6, 5, 4 }, { 16, 14, 12, 10, 8 }, { 24, 21, 18, 15, 12 }, { 32, 28, 24, 20, 16 }, { 40, 35, 30, 25, 20 } };
			double[] testRandomVariable = { 4, 10, 12, 16, 20 };
			double varianceRandomVariable = MatrixCalculations.Covariance(testRandomVariable, testRandomVariable);
			double[,] covarianceMatrix = MatrixCalculations.Covariance(testMatrix); 
			double[,] partialCovTerm = MatrixCalculations.CalculatePartialCovarianceTerm(testMatrix, testRandomVariable, varianceRandomVariable);
			double[,] result = MatrixCalculations.Subtract(covarianceMatrix, partialCovTerm); 
			PrintMatrix(result, true);
		}

		[Test]
		public void TestCalculatePartialCovarianceMatrixAutomatically()
		{
			double[,] testMatrix = { { 8, 7, 6, 5, 4 }, { 16, 14, 12, 10, 8 }, { 24, 21, 18, 15, 12 }, { 32, 28, 24, 20, 16 }, { 40, 35, 30, 25, 20 } };
			double[] testRandomVariable = { 4, 10, 12, 16, 20 };
			double varianceRandomVariable = MatrixCalculations.Covariance(testRandomVariable, testRandomVariable);

			double[,] result = MatrixCalculations.CalculatePartialCovarianceMatrix(testMatrix, testRandomVariable, varianceRandomVariable);
			PrintMatrix(result, true);

		}
	}
}
