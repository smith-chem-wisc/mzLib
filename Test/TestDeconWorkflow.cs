using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using UniDecAPI;
using System.IO;
using MassSpectrometry;
using IO.ThermoRawFileReader;
using System.Runtime.InteropServices;
using System.Diagnostics;
using System.Data;
using OxyPlot; 

namespace Test
{
	unsafe class TestDeconWorkflow
	{
		private MsDataScan scan;
		private Config config;
		private InputUnsafe inp;
		private float[] xarray;
		private float[] yarray;
		public IntPtr isotopeposPtr;
		public IntPtr isotopevalPtr;
		public Decon deconResults;
		public DeconUnsafe deconUnsafe;
		public UnmanagedHandling UnHandler; 

		//[OneTimeSetUp]
		public void Init()
		{
			UnHandler = new UnmanagedHandling(); 
			var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "LowResMS1ScanForDecon.raw");
			//List<MsDataScan> testScan = ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList();
			//scan = testScan[0];

			FilteringParams filterParams = new(numberOfPeaksToKeepPerWindow:1, windowWidthThomsons: 2.006);
			List<MsDataScan> testScanWithFilter = ThermoRawFileReader.LoadAllStaticData(path, filterParams).GetAllScansList();
			scan = testScanWithFilter[0]; 

			// setup the config struct
			UniDecAPIMethods.ConfigMethods.CreateAndSetupConfig(scan, out config);

			// setup the input struct
			inp = UniDecAPIMethods.InputMethods.SetupInputs();

			// assign inp the x and y array data
			xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();			

			inp.dataInt = (float*)UnHandler.AllocateToUnmanaged(xarray);
			inp.dataMZ = (float*)UnHandler.AllocateToUnmanaged(yarray);

			inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
			// correctly assigns a value to inp.barr, but the value isn't correct. Need to make sure that it's 
			// char on the C side. 
			UniDecAPIMethods.UtilityMethods.SetLimits(config, inp);
			IsotopeStruct isoStruct = Isotopes.SetupIsotopes(config, inp);
			config.isolength = Math.Abs(isoStruct.isolength);
			
			inp.isotopeops = (int*)UnHandler.AllocateToUnmanaged(config.isolength * config.lengthmz * config.numz, typeof(int));
			inp.isotopeval = (float*)UnHandler.AllocateToUnmanaged(config.isolength * config.lengthmz * config.numz, typeof(float));

			Isotopes.MakeIsotopesFromAPI(config, inp, isoStruct);
			int numberOfElementsInBarr = config.lengthmz * config.numz;
			deconUnsafe = new DeconUnsafe(); 
			byte[] barr = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.barr, numberOfElementsInBarr);
			/*
			 * inp.barr is a char* pointing to an array of 1s or 0s. These ones and zeroes are stored in an 
			 * 8-bit/1 byte format called ASCII. When converting from C's char to C#'s byte format,
			 * C# is seing the one byte char '49', which is ASCII code for '1' and converting it to 49, instead of 
			 * 1. So we need to convert the ASCII code to the actual number, which means all we need to do to convert
			 * to the correct byte value in C# is subtract the ASCII code for zero, which is 48. The subtraction is implemented in 
			 * the below method. 
			 */
			UniDecAPIMethods.UtilityMethods.ConvertASCIIBytesFromCToByteInCS(ref barr);
			//DirectUniDecPort.Blur.Isotopes.SetupAndMakeIsotopes2(config, inp);

			float threshold = config.psthresh * Math.Abs(config.mzsig) * config.peakshapeinflate;

			int* starttab = (int*) UnHandler.AllocateToUnmanaged(config.lengthmz, typeof(int));
			int* endtab = (int*)UnHandler.AllocateToUnmanaged(config.lengthmz, typeof(int));

			int maxlength = Convolution.SetStartEnds(config, ref inp, starttab, endtab, threshold);

			int pslen = config.lengthmz * maxlength;
			float* mzdist = (float*)UnHandler.AllocateToUnmanaged(pslen, typeof(float));
			float* rmzdist = (float*)UnHandler.AllocateToUnmanaged(pslen, typeof(int));
			int makereverse = 1;
			// makepeakshape2d is very slow (>20s) 
			MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse, starttab, endtab, maxlength);

			int zlength = 1 + 2 * (int)config.zsig;
			int mlength = 1 + 2 * (int)config.msig;
			int* mind = (int*)UnHandler.AllocateToUnmanaged(mlength, typeof(int));
			float* mdist = (float*)UnHandler.AllocateToUnmanaged(mlength, typeof(int)); 
			for (int i = 0; i < mlength; i++)
			{
				mind[i] = i - (mlength - 1) / 2;
				if (config.msig != 0)
				{
					mdist[i] = (float)(Math.Exp(-Math.Pow((i - (zlength - 1) / 2.0), 2)) / (2.0 * config.zsig * config.zsig));
				}
				else
				{
					mdist[i] = 1;
				}
			}
			int* zind = stackalloc int [zlength];
			float* zdist = stackalloc float[zlength];

			int numclose = mlength * zlength;
			int* closemind = (int*)UnHandler.AllocateToUnmanaged(numclose, typeof(int));
			int* closezind = (int*)UnHandler.AllocateToUnmanaged(numclose, typeof(int));
			float* closeval = (float*)UnHandler.AllocateToUnmanaged(numclose, typeof(int));
			int* closeind = (int*)UnHandler.AllocateToUnmanaged(numclose * config.lengthmz * config.numz, typeof(int));
			float* closearray = (float*)UnHandler.AllocateToUnmanaged(numclose * config.lengthmz * config.numz, typeof(int));

			for (int k = 0; k < numclose; k++)
			{
				closemind[k] = mind[k % mlength];
				closezind[k] = zind[(int)k / mlength];
				closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength];
			}

			Normalization.simp_norm_sum(mlength, mdist);
			Normalization.simp_norm_sum(zlength, zdist);
			Normalization.simp_norm_sum(numclose, closeval);

			Blur._MakeSparseBlur(numclose, inp.barr, closezind, closemind, inp.mtab, inp.nztab, 
				inp.dataMZ, closeind, closeval, closearray, config);

			int badness = 1;
			for (int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if (barr[i] == 1)
				{
					badness = 0;
				}
			}
			if (badness == 1)
			{
				throw new InvalidOperationException("Badness = 1...");
			}

			float dmax = MathUtilities.Max(inp.dataInt, config.lengthmz);
			float betafactor = 1;
			if (dmax > 1)
			{
				betafactor = dmax;
			}

			UniDecAPIMethods.UtilityMethods.KillB_CS(inp.dataInt, ref barr, config.intthresh, config.lengthmz, config.numz, config.isolength,
				inp.isotopeops, inp.isotopeval);
			// DirectUniDecPort.FitFunctions.KillB(inp.dataInt, barr, 
			//	config.intthresh, config.lengthmz, config.numz, 
			//	config.isolength, inp.isotopeops, inp.isotopeval); 

			deconUnsafe.blur = (float*)UnHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			deconUnsafe.newblur = (float*)UnHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			float* oldblur = stackalloc float[config.lengthmz * config.numz];
			deconUnsafe.baseline = (float*)UnHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			deconUnsafe.noise = (float*)UnHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float)); 

			for (int i = 0; i < config.lengthmz; i++)
			{
				float val = inp.dataInt[i] / ((float)(config.numz + 2));
				if (config.baselineflag == 1)
				{

					deconUnsafe.baseline[i] = val;
					deconUnsafe.noise[i] = val;
				}

				for (int j = 0; j < config.numz; j++)
				{
					if (barr[ArrayIndexing.Index2D(config.numz, i, j)] == 1)
					{
						if (config.isotopemode == 0)
						{
							deconUnsafe.blur[ArrayIndexing.Index2D(config.numz, i, j)] = val;
						}
						else {deconUnsafe.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 1; }
					}
					else
					{
						deconUnsafe.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 0;
					}
				}
			}

			// copy decon.blur to oldblur: 
			oldblur = deconUnsafe.blur;
			// copy decon.newblur to decon.blur: 
			deconUnsafe.newblur = deconUnsafe.blur;

			float[] dataInt2 = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.dataInt, config.lengthmz);


			Convolution.DeconvolveBaseline(config, inp, deconUnsafe);
			float blurmax = 0F;
			deconUnsafe.conv = 0F;
			
			fixed(float* dataInt2Ptr = &dataInt2[0])
			{
				Blur.PerformIterations(ref deconUnsafe, config, inp, betafactor, maxlength,
					starttab, endtab, mzdist, numclose, closeind, closearray, zlength, mlength,
					closemind, closezind, mdist, dataInt2Ptr, zdist, barr, rmzdist, oldblur);
			}


			if (config.peakshapeinflate != 1 && config.mzsig != 0)
			{
				if (config.speedyflag == 0)
				{
					MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse: 0, starttab, endtab, maxlength);
				}
				else
				{
					MZPeak.MakePeakShape1D(config, inp, threshold, mzdist, rmzdist, makeReverse: 0);
				}
			}

			blurmax = MathUtilities.Max(deconUnsafe.blur, config.lengthmz * config.numz);
			float cutoff = 0F;
			if (blurmax != 0)
			{
				cutoff = 0.000001F;
			}

			ArrayIndexing.ApplyCutoff1D(deconUnsafe.blur, blurmax * cutoff, config.lengthmz * config.numz);

			deconUnsafe.fitdat = (float*)UnHandler.AllocateToUnmanaged(config.lengthmz, typeof(float));

			
			deconUnsafe.error = ErrorFunctions.ErrFunctionSpeedy(config, ref deconUnsafe, barr, inp.dataInt, maxlength, inp.isotopeops, inp.isotopeval,
				starttab, endtab, mzdist); 
			//deconUnsafe.error = ErrorFunctions.ErrFunctionSpeedy(config, deconUnsafe, barr, inp.dataInt, maxlength,
			//	inp.isotopeops, inp.isotopeval, starttab, endtab, mzdist);

			if (config.intthresh != -1)
			{
				for (int i = 0; i < config.lengthmz - 1; i++)
				{
					if (inp.dataInt[i] == 0 && inp.dataInt[i + 1] == 0)
					{
						deconUnsafe.fitdat[i] = 0F;
						deconUnsafe.fitdat[i + 1] = 0F;
					}
				}
			}
			// not tested yet. 
			if (config.isotopemode == 2)
			{
				Isotopes.MonoisotopicToAverageMass(config, inp, deconUnsafe, barr);
			}

			float newblurmax = blurmax;
			if (config.rawflag == 0 || config.rawflag == 2)
			{
				if (config.mzsig != 0)
				{
					newblurmax = Convolution._Reconvolve(config.lengthmz, config.numz, maxlength,
						starttab, endtab, mzdist, deconUnsafe.blur, deconUnsafe.newblur, config.speedyflag, inp.barr);
				}
				else
				{
					deconUnsafe.newblur = deconUnsafe.blur;
				}
			}
			float massmax = config.masslb;
			float massmin = config.massub;

			if (config.fixedmassaxis == 0)
			{
				for (int i = 0; i < config.lengthmz; i++)
				{
					for (int j = 0; j < config.numz; j++)
					{
						if (deconUnsafe.newblur[ArrayIndexing.Index2D(config.numz, i, j)] * barr[ArrayIndexing.Index2D(config.numz, i, j)] > newblurmax * cutoff)
						{
							float testmax = inp.mtab[ArrayIndexing.Index2D(config.numz, i, j)] + threshold * inp.nztab[j] + config.massbins;
							float testmin = inp.mtab[ArrayIndexing.Index2D(config.numz, i, j)] - threshold * inp.nztab[j];

							//To prevent really wierd decimals
							testmin = (float)Math.Round(testmin / config.massbins) * config.massbins;
							testmax = (float)Math.Round(testmax / config.massbins) * config.massbins;

							if (testmax > massmax) { massmax = testmax; }
							if (testmin < massmin) { massmin = testmin; }
						}
					}
				}
			}
			else { massmax = config.massub; massmin = config.masslb; }

			//Checks to make sure the mass axis is good and makes a dummy axis if not
			deconUnsafe.mlen = (int)((massmax - massmin) / config.massbins);
			if (deconUnsafe.mlen < 1)
			{
				massmax = config.massub;
				massmin = config.masslb;
				deconUnsafe.mlen = (int)((massmax - massmin) / config.massbins);

				//Declare the memory

				deconUnsafe.massaxis = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float));
				deconUnsafe.massaxisval = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float));
				deconUnsafe.massgrid = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen * config.numz, typeof(float));

				//Create the mass axis
				for (int i = 0; i < deconUnsafe.mlen; i++)
				{
					deconUnsafe.massaxis[i] = massmin + i * config.massbins;
				}
				deconUnsafe.uniscore = 0;
			}
			else
			{

				//Declare the memory
				deconUnsafe.massaxis = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float));
				deconUnsafe.massaxisval = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float));
				deconUnsafe.massgrid = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen * config.numz, typeof(float));


				//Create the mass axis
				for (int i = 0; i < deconUnsafe.mlen; i++)
				{
					deconUnsafe.massaxis[i] = massmin + i * config.massbins;
				}
			}

			MassIntensityDetermination.IntegrateMassIntensities(config, ref deconUnsafe, inp);
		}
		//[OneTimeTearDown]
		public void TearDown()
		{
			foreach(IntPtr i in UnHandler.ListPtrs)
			{
				Marshal.FreeHGlobal(i);
			}
		}

		[Test]
		public void TestUniDecDeconvolutionWorklow()
		{
			float scorethresh = 0f;
			config.peakwin = 20;

			deconUnsafe.peakx = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float));
			deconUnsafe.peaky = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float)); 
			
			deconUnsafe.uniscore = Scoring.UniScorePorted(config, ref deconUnsafe, inp, scorethresh, config.peakwin, UnHandler);
			PrintProperties(deconUnsafe); 
		}
		[Test]
		public void TestDeconResults()
		{
			float scorethresh = 0f;
			config.peakwin = 20;

			deconUnsafe.peakx = (float*)UnHandler.AllocateToUnmanaged(deconResults.mlen, typeof(float));
			deconUnsafe.peaky = (float*)UnHandler.AllocateToUnmanaged(deconResults.mlen, typeof(float));

			deconUnsafe.uniscore = Scoring.UniScorePorted(config, ref deconUnsafe, inp, scorethresh, config.peakwin, UnHandler);
			DeconResults deconFin = new(deconUnsafe, config);
			Console.WriteLine(string.Join(", ", deconFin.PeakX.ToList())); 
		}

		[Test]
		[TestCase(1000.3f)]
		[TestCase(8129f)]
		[TestCase(9087.12f)]
		[TestCase(1.12f)]
		[TestCase(85.0f)]
		public void TestNearFastCS(float searchVal)
		{

			int numData = 100000; 
			float[] seq = Enumerable.Range(1, numData).Select(i => i * 1f).ToArray();
			// returns the index of the closest point.
			// int nearFastResult = DirectUniDecPort.Convolution.NearFastCS(seq, 50f, numData);
			// Assert.AreEqual(49, nearFastResult); 

			// total time taken for 100 points: 35 ms. So it's actually extremely slow to be searching all these
			// points with the algorithm as-is. 

			// binary search recursive: 
			Stopwatch stopwatch = new();
			
			stopwatch.Start();
			int binarySearchResultIndex = BinarySearch(seq, searchVal, 0, numData - 1);
			stopwatch.Stop();

			Assert.AreEqual((int)searchVal - 1, binarySearchResultIndex);
			Console.WriteLine("binary search with recursion: " + stopwatch.ElapsedTicks.ToString());
			
			fixed(float* seqPtr = &seq[0])
			{
				Stopwatch stpwtch = new();
				stpwtch.Start(); 
				int resultNearFast = Convolution.Nearfast(seqPtr, searchVal, numData);
				stpwtch.Stop();
				Assert.AreEqual((int)searchVal - 1, resultNearFast); 
				Console.WriteLine("nearfast search: " + stpwtch.ElapsedTicks.ToString()); 
			}

		}
		[Test]
		public void TestSetupIsotopes()
		{
			IsotopeStruct isoStruct = Isotopes.SetupIsotopes(config, inp);
			PrintProperties(isoStruct); 
			for(int i = 0; i < config.lengthmz; i++)
			{
				for(int j = 0; j < config.numz; j++)
				{
					for(int k = 0; k < 1; k++)
					{
						Console.WriteLine(inp.isotopeval[ArrayIndexing.Index3D(config.numz, 1, i, j, k)].ToString()); 
					}
				}
			}
		}
		[Test]
		public void TestMakeIsotopes()
		{
			IsotopeStruct isoStruct = Isotopes.SetupIsotopes(config, inp);
			config.isolength = Math.Abs(isoStruct.isolength); 
			isotopeposPtr = Marshal.AllocHGlobal(Marshal.SizeOf(1) * config.isolength * config.lengthmz * config.numz);
			isotopevalPtr = Marshal.AllocHGlobal(Marshal.SizeOf(1f) * config.isolength * config.lengthmz * config.numz);
			inp.isotopeops = (int*)isotopeposPtr;
			inp.isotopeval = (float*)isotopevalPtr;
			Isotopes.MakeIsotopesFromAPI(config, inp, isoStruct);
		}
		[Test]
		public void TestIndex3DArray()
		{
			int[,,] test3DArray = new int[,,] 
			{
				{ { 0, 1, 2 } },
				{ { 3, 4, 5 } },
				{ { 6, 7, 8 } } 
			};
			int[] test1DArray = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
			int result = ArrayIndexing.Index3D(1, 3, 0, 1, 2);
			int result2 = ArrayIndexing.Index3D(1, 3, 1, 2, 0); 
			Assert.AreEqual(test1DArray[5], test1DArray[result]);
			Assert.AreEqual(test1DArray[7], test1DArray[result2]);
			//int shouldFail = test1DArray[DirectUniDecPort.ArrayIndexing.Index3D(1, 3, 3, 3, 3)]; 
		}

		public static int BinarySearch(float[] array, float searchVal, int leftIndex, int rightIndex)
		{
			if (rightIndex >= leftIndex)
			{
				int mid = leftIndex + (rightIndex - leftIndex) / 2;

				if (array[mid] == searchVal) return mid;
				if (array[mid] > searchVal) return BinarySearch(array, searchVal, leftIndex, mid - 1);
				return BinarySearch(array, searchVal, mid + 1, rightIndex); 
			}
			if (Math.Abs(searchVal - array[leftIndex]) >= Math.Abs(searchVal - array[rightIndex]))
			{
				return rightIndex;
			}
			else
			{
				return leftIndex; 
			}
		}
		public void PrintProperties(object o)
		{
			foreach (var field in o.GetType().GetFields())
			{
				Console.WriteLine(field.Name + ": " + field.GetValue(o));
			}
		}
		public void TestBinning()
		{
			
		}
		[Test]
		[TestCase(0.1f)]
		[TestCase(0.01f)]
		[TestCase(0f)]
		public void TestPeakDetect(float threshold)
		{
			deconResults.peakx = new float[deconResults.mlen];
			deconResults.peaky = new float[deconResults.mlen];
			int result = 0; 
			fixed (float* peakxPtr = &deconResults.peakx[0], peakyPtr = &deconResults.peaky[0])
			{
				result = Scoring.PeakDetect(inp.dataMZ, inp.dataInt, config.lengthmz, 20,
					threshold, peakxPtr, peakyPtr);
			}
			Console.WriteLine(result.ToString()); 
		}
		[Test]
		public void TestGetFWHMS()
		{
			float threshold = 0.01f;
			deconUnsafe.peakx = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float));
			deconUnsafe.peaky = (float*)UnHandler.AllocateToUnmanaged(deconUnsafe.mlen, typeof(float)); 
			// get the number of peaks: 
			int plen = Scoring.PeakDetect(inp.dataMZ, inp.dataInt, config.lengthmz, 20,
				threshold, deconUnsafe.peakx, deconUnsafe.peaky);
			float* fwhmHigh = stackalloc float[plen];
			float* fwhmLow = stackalloc float[plen];
			float* badFwhm = stackalloc float[plen]; 

			Scoring.GetFWHMS(config, plen, deconResults.mlen, deconUnsafe.massaxis, deconUnsafe.massaxisval, deconUnsafe.peakx,
				fwhmLow, fwhmHigh, badFwhm);
			Console.WriteLine(string.Join("; ", fwhmLow[0], fwhmHigh[0], badFwhm[0]));
			
		}
		[Test]
		[TestCase(1, 2.006)]
		[TestCase(5, 2.006)]
		[TestCase(1, 1.003)]
		[TestCase(5, 1.003)]

		public void TestUniDecDeconEngine(int peakPerWindow, double windowWidth)
		{
			var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "LowResMS1ScanForDecon.raw");
			List<MsDataScan> testScan = ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList();
			scan = testScan[0];

			FilteringParams filterParams = new(numberOfPeaksToKeepPerWindow: peakPerWindow, windowWidthThomsons: windowWidth);
			List<MsDataScan> testScanWithFilter = ThermoRawFileReader.LoadAllStaticData(path, filterParams).GetAllScansList();
			scan = testScanWithFilter[0];
			
			scan.MassSpectrum.DeconvoluteWithUniDec(out DeconResults deconResults);
			Console.WriteLine(string.Join("\n", deconResults.DScores.ToList())); 
		}
		[Test]
		public void TestGenerateGaussianKernel()
		{
			UniDecDeconEngine engine = new();
			double[] result = engine.GenerateGuassianKernel(1, 5);
			double[] expected = new double[] { 0.053991, 0.2419707, 0.39894, 0.24197, 0.053991};
			for(int i = 0; i < result.Length; i++)
			{
				Assert.AreEqual(expected[i], result[i], 0.001); 
			}
		}
		[Test]
		public void TestConvoluteFull()
		{
			UniDecDeconEngine engine = new(); 

			double[] testSignal = new double[] { 0, 0, 0, 0, 255, 255, 255, 0, 0, 0, 0 };
			double[] kernel = engine.GenerateGuassianKernel(0.5, 3);
			
			double[] result = engine.ConvoluteFull(testSignal, kernel);
			Console.WriteLine(string.Join("\n", result.ToList())); 
			// TODO: Write and actual assert statement here
		}

		[Test]
		public void TestCreateChargeAndMZMatrix()
		{
			UniDecDeconEngine engine = new(); 
			string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "UniDecDeconTestData.txt");
			DataTable proteinDt = LoadTsvFileToDataTable(path);

			var testMzData = proteinDt.AsEnumerable().Select(x => new
			{
				Mz = (double)x["Mz"], 
				Int = (double)x["Intensity"]
			}).ToList();

			double[] xarray = testMzData.Select(i => (double)i.Mz).ToArray();
			double[] yarray = testMzData.Select(i => (double)i.Int).ToArray();

			double[,] matrix = engine.CreateChargeAndMZMatrix(yarray, 5, 50);

			int colLength = matrix.GetLength(1); 
			for(int i = 0; i < matrix.GetLength(0); i++)
			{
				Console.WriteLine(string.Join(", ", Enumerable.Range(0, matrix.GetLength(1))
					.Select(x => matrix[i, x]).ToArray())); 
			}
		}
		[Test]
		public void TestDeconvolutionII()
		{
			CreateBetterTestData(out double[] yarray, out double[] xarray); 
			// setup
			UniDecDeconEngine engine = new();

			double[,] measuredSpectrum = new double[xarray.Length, 2]; 
			for(int i = 0; i < xarray.Length; i++)
			{
				measuredSpectrum[i, 0] = xarray[i];
				measuredSpectrum[i, 1] = yarray[i]; 
			}

			/*
			 * iteration = delta functions * measured spectrum / charge axis convolved with peak shape 
			 * 
			 * delta functions are smoothed by convoluting with a filter before proceeding with the deconvolution
			 **/

			// creates a deep copy of matrix and creates the result of the iterations. 
			double[,] result = engine.CreateDeepCopy(measuredSpectrum);
			// create the gaussian kernel. 
			// spacing for this example is 0.1 between each measured m/z. 
			// so the minimum standard deviation of a kernel is going to be 2 * 0.1. 
			// This makes sense since you can't meaningfully convolute a signal with a function that is the min frequency/2. 
			double[] gaussKernel1D = engine.GenerateGuassianKernel(0.2, 8);

			
			// f^(t=0) = intensity data 
			double[,] initialFt = engine.CreateChargeAndMZMatrix(yarray, 5, 50);
			initialFt = engine.ApplyLogMeanFilter(initialFt, 7); 
			double[] initialSummedChargeAxis = engine.SumChargeAxis(initialFt);
			double[] initialConvolutedChargeAxis = engine.ConvoluteSummedChargeAxis(initialSummedChargeAxis, gaussKernel1D);
			double[] invConvChargeAxis = engine.TakeReciprocalOfArray(initialConvolutedChargeAxis);

			double[] temp = engine.ElementwiseMultiplyArrays(yarray, invConvChargeAxis);
			double[] iterationResult = engine.ElementwiseMultiplyArrays(yarray, temp);
			
			double[,] ft = new double[45, 10000];
			double[,] smoothedIterationResult = new double[45, 10000]; 
			double[] summedChargeAxis = new double[10000];
			double[] convolutedChargeAxis = new double[10000];
			double[] recipChargeAxis = new double[10000]; 
			
			for (int i = 0; i < 10; i++)
			{
				ft = engine.CreateChargeAndMZMatrix(iterationResult, 5, 50); 
				smoothedIterationResult = engine.ApplyLogMeanFilter(ft, 7);
				summedChargeAxis = engine.SumChargeAxis(smoothedIterationResult);
				convolutedChargeAxis = engine.ConvoluteSummedChargeAxis(summedChargeAxis, gaussKernel1D);
				// h(x) / c(x)
				// reciprocal followed by multiplication should hypothetically be faster than element-wise division of the two arrays. 
				recipChargeAxis = engine.TakeReciprocalOfArray(convolutedChargeAxis);
				temp = engine.ElementwiseMultiplyArrays(yarray, recipChargeAxis);

				// f^t * h(x) / c(x)
				iterationResult = engine.ElementwiseMultiplyArrays(iterationResult, temp); 
			}

			int[] chargeArray = Enumerable.Range(5, 45).ToArray();
			double[,] massTable = engine.CreateChargeByMzTable(chargeArray, xarray, 1.007);

			engine.CreateMassAxes(out double[] massaxis, out double[] massAxisVals, 50000, 10000, 5);
			double[,] deconMassTable = engine.CreateDeconvolutedMassTable(massaxis, chargeArray); 
			engine.IntegrateTransform(massTable, massaxis, ref massAxisVals, smoothedIterationResult, ref deconMassTable);

			// PrintMatrix(deconMassTable);
			Console.WriteLine(string.Join("\n", massAxisVals.AsEnumerable())); 
			double[] deconMassSpectrum = engine.ColSums(deconMassTable);
			CreateDeconMassSpectrum(deconMassSpectrum, massaxis, "UniDecOuptut.pdf"); 
		}
		public void CreateDeconMassSpectrum(double[] deconMassSpec, double[] massAxis, string fileName)
		{
			var data = new OxyPlot.Series.LineSeries()
			{
				Title = $"Deconvoluted Mass Spectra",
				Color = OxyPlot.OxyColors.Blue,
				StrokeThickness = 1,
				MarkerSize = 1,
				MarkerType = OxyPlot.MarkerType.Circle
			}; 
			for(int i = 0; i < deconMassSpec.Length; i++)
			{
				data.Points.Add(new OxyPlot.DataPoint(massAxis[i], deconMassSpec[i])); 
			}
			var model = new OxyPlot.PlotModel
			{
				Title = $"Deconvoluted Mass Spectrum"
			};

			model.Series.Add(data);
			WriteOxyPlotToPDF(fileName, model); 
		}
		public static void WriteOxyPlotToPDF(string fileName, PlotModel plotmodel)
		{
			using (var stream = File.Create(fileName))
			{
				var pdfExporter = new PdfExporter { Width = 600, Height = 600 };
				pdfExporter.Export(plotmodel, stream);
			}
		}
		public void PrintMatrix(double[,] matrix)
		{
			for(int i = 0; i < matrix.GetLength(0); i++)
			{
				Console.WriteLine(string.Join(", ", Enumerable.Range(0, matrix.GetLength(1))
					.Select(x => matrix[i, x]).ToArray())); 
			}
		}
		[Test]
		public void TestCreateMassAxes()
		{
			CreateBetterTestData(out double[] yarray, out double[] xarray);
			double massMax = 50000;
			double massMin = 10000;
			double massBinWidth = 5;

			UniDecDeconEngine eng = new();
			eng.CreateMassAxes(out double[] massAxis, out double[] massAxisVals, massMax, massMin, massBinWidth);

			Console.WriteLine(string.Join("\n", massAxis.AsEnumerable())); 

		}
		[Test]
		public void TestCreateChargeByMzTable()
		{
			UniDecDeconEngine eng = new();
			CreateBetterTestData(out double[] yarray, out double[] xarray);
			int[] chargeArray = Enumerable.Range(5, 45).ToArray();
			double[,] massTable = eng.CreateChargeByMzTable(chargeArray, xarray, 1.007); 
		}
		public DataTable LoadTsvFileToDataTable(string tsvPath)
		{
			DataTable dt = new DataTable();
			using (StreamReader reader = new StreamReader(tsvPath))
			{
				string line = "";
				int rowcount = 0;
				while ((line = reader.ReadLine()) != null)
				{
					string[] tsv = line.Split("\t").ToArray();
					if (rowcount == 0)
					{
						foreach (string colName in tsv)
						{
							dt.Columns.Add(colName, typeof(double));
						}
						rowcount++; 
						continue;
					}

					dt.Rows.Add(tsv);
					rowcount++; 
				}
			}
			return dt; 
		}
		public void CreateBetterTestData(out double[] intArray, out double[] mzArray)
		{
			double highMz = 1500;
			double lowMz = 500;
			double deltaMz = 0.1;
			int totalValues = (int)((highMz - lowMz) / deltaMz);
			double[] mzAxis = Enumerable.Range(0, totalValues).Select(i => i*deltaMz + lowMz).ToArray();
			double[] intensityArray = new double[totalValues]; 
			List<double> mzToHaveIntensity = (new double[] { 1201.0, 1091.9, 1001.0, 924.1, 858.2 }).ToList();
			mzToHaveIntensity.Sort(); 

			int[] indexes = new int[mzToHaveIntensity.Count];

			for(int i = 0; i < mzToHaveIntensity.Count; i++)
			{
				indexes[i] = Array.BinarySearch(mzAxis, mzToHaveIntensity[i]);
				if(indexes[i] < 0)
				{
					indexes[i] = ~indexes[i];
				}
			}
			
			var indexesList = indexes.ToList();
			int j = 0;
			double[] intensityValues = new double[] { 0.0309, 0.0619, 0.0929, 0.124, 0.0929 }; 
			
			for(int i = 0; i < intensityArray.Length; i++)
			{
				intensityArray[i] = 0; 

				if (indexes[j] == i)
				{
					intensityArray[i] = intensityValues[j]; 
					j++;
					// loop will try to evaluate beyond indexes.Length, so 
					// needed to break the loop once last indexes[j] is reached. Also saves time. 
					if(j == indexes.Length)
					{
						break; 
					}
				}
			}
			intArray = intensityArray;
			mzArray = mzAxis; 
		}
	}
}

