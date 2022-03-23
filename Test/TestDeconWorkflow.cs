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

		[OneTimeSetUp]
		public void Init()
		{
			var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "LowResMS1ScanForDecon.raw");
			List<MsDataScan> testScan = ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList();
			scan = testScan[0];

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
		}
		[OneTimeTearDown]
		public void TearDown()
		{
			Marshal.FreeHGlobal(xarrayPtr);
			Marshal.FreeHGlobal(yarrayPtr);
		}

		[Test]
		public void TestUniDecDeconvolutionWorklow()
		{
			int numberOfElementsInBarr = config.lengthmz * config.numz;
			
			char[] barr = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.barr, numberOfElementsInBarr);
			float threshold = config.psthresh * Math.Abs(config.mzsig) * config.peakshapeinflate;

			int[] starttab = new int[config.lengthmz];
			int[] endtab = new int[config.lengthmz];
			int maxlength = DirectUniDecPort.Convolution.SetStartEnds(config, ref inp, ref starttab, ref endtab, threshold);
			
			int pslen = config.lengthmz * maxlength;
			float[] mzdist = new float[pslen];
			float[] rmzdist = new float[pslen];
			int makereverse = 1;
			// makepeakshape2d is very slow (>20s) 
			DirectUniDecPort.MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse, starttab, endtab, maxlength);

			int zlength = 1 + 2 * (int)config.zsig;
			int mlength = 1 + 2 * (int)config.msig;
			int[] mind = new int[mlength];
			float[] mdist = new float[mlength]; 

			for(int i = 0; i < mlength; i++)
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

			for(int k = 0; k < numclose; k++)
			{
				closemind[k] = mind[k % mlength];
				closezind[k] = zind[(int)k / mlength];
				closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength]; 
			}

			DirectUniDecPort.Normalization.SimpNormSum(mlength, mdist);
			DirectUniDecPort.Normalization.SimpNormSum(zlength, zdist);
			DirectUniDecPort.Normalization.SimpNormSum(numclose, closeval);

			DirectUniDecPort.Blur.MakeSparseBlur(inp, config, numclose, barr, closezind,
				closemind, closeind, closeval, closearray);

			int badness = 1; 
			for(int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if(barr[i] == 1)
				{
					badness = 0; 
				}
			}
			if(badness == 1)
			{
				throw new InvalidOperationException("Badness = 1..."); 
			}

			// return Decon 
		}
	}
}
