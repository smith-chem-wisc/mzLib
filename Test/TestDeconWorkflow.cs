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

namespace Test
{
	unsafe class TestDeconWorkflow
	{
		private MsDataScan scan;
		private Config config;
		private InputUnsafe inp;
		private float[] xarray;
		private float[] yarray;
		public IntPtr xarrayPtr;
		public IntPtr yarrayPtr;
		public IntPtr isotopeposPtr;
		public IntPtr isotopevalPtr;
		public Decon deconResults;
		public DeconUnsafe deconUnsafe; 

		[OneTimeSetUp]
		public void Init()
		{
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



			xarrayPtr = Marshal.AllocHGlobal(Marshal.SizeOf(xarray[0]) * xarray.Length);
			Marshal.Copy(xarray, 0, xarrayPtr, xarray.Length); 

			yarrayPtr = Marshal.AllocHGlobal(Marshal.SizeOf(yarray[0]) * yarray.Length);
			Marshal.Copy(yarray, 0, yarrayPtr, yarray.Length); 

			

			inp.dataInt = (float*)yarrayPtr;
			inp.dataMZ = (float*)xarrayPtr;

			inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
			// correctly assigns a value to inp.barr, but the value isn't correct. Need to make sure that it's 
			// char on the C side. 
			UniDecAPIMethods.UtilityMethods.SetLimits(config, inp);
			IsotopeStruct isoStruct = Isotopes.SetupIsotopes(config, inp);
			config.isolength = Math.Abs(isoStruct.isolength);
			
			isotopeposPtr = Marshal.AllocHGlobal(Marshal.SizeOf(1) * config.isolength * config.lengthmz * config.numz);
			isotopevalPtr = Marshal.AllocHGlobal(Marshal.SizeOf(1f) * config.isolength * config.lengthmz * config.numz);
			inp.isotopeops = (int*)isotopeposPtr;
			inp.isotopeval = (float*)isotopevalPtr;
			Isotopes.MakeIsotopesFromAPI(config, inp, isoStruct);
			int numberOfElementsInBarr = config.lengthmz * config.numz;
			deconResults = new Decon();
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

			int[] starttab = new int[config.lengthmz];
			int[] endtab = new int[config.lengthmz];
			int maxlength = Convolution.SetStartEnds(config, ref inp, ref starttab, ref endtab, threshold);

			int pslen = config.lengthmz * maxlength;
			float[] mzdist = new float[pslen];
			float[] rmzdist = new float[pslen];
			int makereverse = 1;
			// makepeakshape2d is very slow (>20s) 
			MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse, starttab, endtab, maxlength);

			int zlength = 1 + 2 * (int)config.zsig;
			int mlength = 1 + 2 * (int)config.msig;
			int[] mind = new int[mlength];
			float[] mdist = new float[mlength];

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
			int[] zind = new int[zlength];
			float[] zdist = new float[zlength];

			int numclose = mlength * zlength;
			int[] closemind = new int[numclose];
			int[] closezind = new int[numclose];
			float[] closeval = new float[numclose];
			int[] closeind = new int[numclose * config.lengthmz * config.numz];
			float[] closearray = new float[numclose * config.lengthmz * config.numz];

			for (int k = 0; k < numclose; k++)
			{
				closemind[k] = mind[k % mlength];
				closezind[k] = zind[(int)k / mlength];
				closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength];
			}

			Normalization.SimpNormSum(mlength, mdist);
			Normalization.SimpNormSum(zlength, zdist);
			Normalization.SimpNormSum(numclose, closeval);

			Blur.MakeSparseBlur(inp, config, numclose, barr, closezind,
				closemind, closeind, closeval, closearray);

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

			deconUnsafe.blur = (float*)UnmanagedHandling.AllocateToUnmanaged(config.lengthmz * config.numz);
			deconUnsafe.newblur = (float*)UnmanagedHandling.AllocateToUnmanaged(config.lengthmz * config.numz);
			float* oldblur = (float*)UnmanagedHandling.AllocateToUnmanaged(config.lengthmz * config.numz);
			deconUnsafe.baseline = (float*)UnmanagedHandling.AllocateToUnmanaged(config.lengthmz * config.numz);
			deconUnsafe.noise = (float*)UnmanagedHandling.AllocateToUnmanaged(config.lengthmz * config.numz); 

			deconResults.blur = new float[config.lengthmz * config.numz];
			deconResults.newblur = new float[config.lengthmz * config.numz];
			//float[] oldblur = new float[config.lengthmz * config.numz];
			deconResults.baseline = new float[config.lengthmz * config.numz];
			deconResults.noise = new float[config.lengthmz * config.numz];

			for (int i = 0; i < config.lengthmz; i++)
			{
				float val = inp.dataInt[i] / ((float)(config.numz + 2));
				if (config.baselineflag == 1)
				{

					deconResults.baseline[i] = val;
					deconResults.noise[i] = val;
				}

				for (int j = 0; j < config.numz; j++)
				{
					if (barr[ArrayIndexing.Index2D(config.numz, i, j)] == 1)
					{
						if (config.isotopemode == 0)
						{
							deconResults.blur[ArrayIndexing.Index2D(config.numz, i, j)] = val;
						}
						else { deconResults.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 1; }
					}
					else
					{
						deconResults.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 0;
					}
				}
			}

			// copy decon.blur to oldblur: 
			oldblur = deconResults.blur;
			// copy decon.newblur to decon.blur: 
			deconResults.newblur = deconResults.blur;

			float[] dataInt2 = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.dataInt, config.lengthmz);

			Convolution.DeconvolveBaseline(config, inp, deconResults);
			float blurmax = 0F;
			deconResults.conv = 0F;
			int off = 0;

			Blur.PerformIterations(ref deconResults, config, inp, betafactor, maxlength,
				starttab, endtab, mzdist, numclose, closeind, closearray, zlength, mlength,
				closemind, closezind, mdist, dataInt2, zdist, barr, rmzdist, oldblur);


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

			blurmax = MathUtilities.Max(deconResults.blur, config.lengthmz * config.numz);
			float cutoff = 0F;
			if (blurmax != 0)
			{
				cutoff = 0.000001F;
			}

			ArrayIndexing.ApplyCutoff1D(ref deconResults.blur, blurmax * cutoff, config.lengthmz * config.numz);

			deconResults.fitdat = new float[config.lengthmz];

			// This line of code is currently not working. 
			deconResults.error = ErrorFunctions.ErrFunctionSpeedy(config, deconResults, barr, inp.dataInt, maxlength,
				inp.isotopeops, inp.isotopeval, starttab, endtab, mzdist);

			if (config.intthresh != -1)
			{
				for (int i = 0; i < config.lengthmz - 1; i++)
				{
					if (inp.dataInt[i] == 0 && inp.dataInt[i + 1] == 0)
					{
						deconResults.fitdat[i] = 0F;
						deconResults.fitdat[i + 1] = 0F;
					}
				}
			}
			// not tested yet. 
			if (config.isotopemode == 2)
			{
				Isotopes.MonoisotopicToAverageMass(config, inp, deconResults, barr);
			}

			float newblurmax = blurmax;
			if (config.rawflag == 0 || config.rawflag == 2)
			{
				if (config.mzsig != 0)
				{
					newblurmax = Convolution.Reconvolve(config.lengthmz, config.numz, maxlength,
						starttab, endtab, mzdist, deconResults.blur, deconResults.newblur, config.speedyflag, barr);
				}
				else
				{
					deconResults.newblur = deconResults.blur;
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
						if (deconResults.newblur[ArrayIndexing.Index2D(config.numz, i, j)] * barr[ArrayIndexing.Index2D(config.numz, i, j)] > newblurmax * cutoff)
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
			deconResults.mlen = (int)((massmax - massmin) / config.massbins);
			if (deconResults.mlen < 1)
			{
				massmax = config.massub;
				massmin = config.masslb;
				deconResults.mlen = (int)((massmax - massmin) / config.massbins);

				//Declare the memory

				deconResults.massaxis = new float[deconResults.mlen];
				deconResults.massaxisval = new float[deconResults.mlen];
				deconResults.massgrid = new float[deconResults.mlen * config.numz];

				//Create the mass axis
				for (int i = 0; i < deconResults.mlen; i++)
				{
					deconResults.massaxis[i] = massmin + i * config.massbins;
				}
				deconResults.uniscore = 0;
			}
			else
			{

				//Declare the memory
				deconResults.massaxis = new float[deconResults.massaxis.Length];
				deconResults.massaxisval = new float[deconResults.mlen];
				deconResults.massgrid = new float[deconResults.mlen * config.numz];


				//Create the mass axis
				for (int i = 0; i < deconResults.mlen; i++)
				{
					deconResults.massaxis[i] = massmin + i * config.massbins;
				}
			}

			MassIntensityDetermination.IntegrateMassIntensities(config, deconResults, inp);
		}
		[OneTimeTearDown]
		public void TearDown()
		{
			Marshal.FreeHGlobal(xarrayPtr);
			Marshal.FreeHGlobal(yarrayPtr);
			Marshal.FreeHGlobal(isotopevalPtr);
			Marshal.FreeHGlobal(isotopeposPtr); 
		}

		[Test]
		public void TestUniDecDeconvolutionWorklow()
		{
			float scorethresh = 0f;
			config.peakwin = 20; 
			//deconResults.uniscore = Scoring.UniScorePorted(config, deconResults, inp, scorethresh);
			Console.WriteLine(deconResults.peakx[1].ToString() + deconResults.peaky[1].ToString() + deconResults.dscores[1].ToString()); 
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
			fixed (float* massaxisPtr = &deconResults.massaxis[0], massaxisvalPtr = &deconResults.massaxisval[0],
				peakxPtr = &deconResults.peakx[0], peakyPtr = &deconResults.peaky[0], newblurPtr = &deconResults.newblur[0], 
				massgridPtr = &deconResults.massgrid[0])
			{
				// get the number of peaks: 
				int plen = Scoring.PeakDetect(inp.dataMZ, inp.dataInt, config.lengthmz, 20,
					threshold, peakxPtr, peakyPtr);
				float[] fwhmHigh;
				float[] fwhmLow;
				float[] badFwhm; 

				Scoring.GetFWHMS(config, plen, deconResults.mlen, massaxisPtr, massaxisvalPtr, peakxPtr,
					out fwhmLow, out fwhmHigh, out badFwhm);

				Console.WriteLine(string.Join("; ", fwhmLow.Length, fwhmHigh.Length, badFwhm.Length));
				Console.WriteLine(string.Join("; ", fwhmLow[0], fwhmHigh[0], badFwhm[0]));
			}
		}
		[Test]
		public void TestScoreFromPeaksPorted()
		{

		}
	}
}
