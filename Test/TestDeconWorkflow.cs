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
			DirectUniDecPort.MZPeak.MakePeakShape2D(config.lengthmz, maxlength, ref starttab, ref endtab, inp.dataMZ,
				Math.Abs(config.mzsig) * config.peakshapeinflate, config.psfun, config.speedyflag, ref mzdist, ref rmzdist,
				makereverse); 


			// return Decon 
		}
	}
}
