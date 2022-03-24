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
	// minimize your unsafe code calling!!! 
	// If you have to use it, minimize the unsafe scope!! 

	class TestUniDecAPI
	{

		private MsDataScan scan;
		[OneTimeSetUp]
		public void Init()
		{
			var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "LowResMS1ScanForDecon.raw");
			List<MsDataScan> testScan = ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList();
			scan = testScan[0];
		}
		[Test]
		public void TestModifyConfigCS()
		{
			Config testConfig = new();
			testConfig.adductmass = 1.003F;
			testConfig = UniDecAPIMethods.ConfigMethods.ReturnModifiedConfigToCS(testConfig);
			PrintProperties(testConfig);
		}
		[Test]
		public void TestModifyConfigToDefault()
		{
			// this works. requires some unsafe code. 
			Config testConfig = new();
			testConfig = UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(testConfig);
			PrintProperties(testConfig);
		}
		[Test]
		public void TestDecon()
		{
			Decon decon = new();
			try
			{
				// creates Decon object in C dll code and passes it back to C# 
				decon = UniDecAPIMethods.DeconMethods.SetupDecon();
				// prints the fields of the Decon object in C#. 
				PrintProperties(decon);
			}
			finally
			{
				// frees the Decon object memory from C dll. 
				UniDecAPIMethods.DeconMethods.FreeDecon(decon);
			}
		}
		[Test]
		public void TestInput()
		{
			float[] xarray = { 1.1F, 2.2F, 3.3F };
			float[] yarray = { 2.2F, 3.3F, 4.4F };
			unsafe
			{
				InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();
				fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
				{


					inp.dataMZ = ptrXarray;
					inp.dataInt = ptrYarray;
					Console.WriteLine("Inside 'Fixed' statement: " + inp.dataMZ[0].ToString());
					Assert.AreEqual(2.2F, inp.dataMZ[1]);
				}
				Console.WriteLine("Outside 'Fixed' Statement: " + inp.dataMZ[0].ToString());
			}
		}
		[Test]
		public void TestAverage()
		{
			// failing in unidec minimal
			float[] testArray = { 1.1F, 1.2F, 5.5F, 6.6F };
			float average = UniDecAPIMethods.DeconMethods.Average(testArray.Length, testArray);
			Console.WriteLine(average.ToString());
		}
		public void PrintProperties(object o)
		{
			foreach (var field in o.GetType().GetFields())
			{
				Console.WriteLine(field.Name + ": " + field.GetValue(o));
			}
		}
		[Test]
		public void TestAverageInOtherDll()
		{
			// passing in the testDLL project. 
			// iterate through and find which files are not working. 
			float[] testArray = { 1.1F, 2.2F, 33.3F };
			float result = APITesting.Average(testArray);
			Console.WriteLine(result.ToString());
		}
		[Test]
		public void TestExecuteTestFFT()
		{
			// based on this, the dependency issue is not in the fftw3 library. Though I do need to change 
			// the address to the new one. 
			int result = APITesting.ExecuteTestFFT();
			Assert.AreEqual(0, result);
		}
		[Test]
		public void TestArrayIndexingFindMax()
		{
			float[] array1 = { 1.1F, 100.0F, 121.123F };
			float result = APITesting.Max1d(array1);
			Assert.AreEqual(121.123F, result);
		}
		[Test]
		public void TestFloatIntCase()
		{
			int testInt = 23;
			float castIntToFloat = (float)testInt;
			Console.WriteLine(castIntToFloat);
		}
		[Test]
		public void TestSetupInputs()
		{
			InputUnsafe inp = new();
			try
			{
				inp = UniDecAPIMethods.InputMethods.SetupInputs();
				PrintProperties(inp);
			}
			finally
			{
				UniDecAPIMethods.InputMethods.FreeInputs(inp);
			}
		}

		[Test]
		public void TestSetupConfig()
		{

			float minOfScan = (float)scan.MassSpectrum.XArray.Min();

			Config config = new();
			UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(config);
			UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);
			Assert.AreEqual(minOfScan, config.minmz);
		}
		[Test]
		public unsafe void TestConfigAndInputsSetup()
		{
			Config config = new();
			InputUnsafe inp = new();
			Decon decon = new();

			config = UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(config);
			UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);
			inp = UniDecAPIMethods.InputMethods.SetupInputs();

			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();

			fixed (float* xarrayPtr = &xarray[0], yarrayPtr = &yarray[0])
			{
				inp.dataMZ = xarrayPtr;
				inp.dataInt = yarrayPtr;

				inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
				PrintProperties(inp);
				Console.WriteLine("Now printing the config object" + "\n");
				PrintProperties(config);
				for (int i = 0; i < 100; i++) { Console.WriteLine(inp.barr[i].ToString()); }
			}
		}
		[Test]
		public void TestMemoryAllocationOfBarr()
		{
			int result = AlgorithmTesting.MemoryAllocationOfBarr();
			Assert.AreEqual(1, result);
		}
		[Test]
		public void TestAllocateMemoryToPointersThenFree()
		{
			int result = AlgorithmTesting.AllocateMemoryToPointersThenFree();
			Assert.AreEqual(1, result);
		}
		[Test]
		public void TestUseMemcpyInC()
		{
			char[] expected = "1".ToCharArray();
			char result = AlgorithmTesting.UseMemcpyInC();
			Assert.AreEqual(expected[0], result);
		}
		[Test]
		public void TestUseMemcpyWithInpAndConfig()
		{
			char[] expected = "1".ToCharArray();
			char result = AlgorithmTesting.UseMemcpyWithInpAndConfigObjects();
			Assert.AreEqual(expected[0], result);
		}
		[Test]
		public void TestCFuncTestSEtStartEnds()
		{
			Config config = new();
			InputUnsafe inp = new();

			config = UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(config);
			UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);
			inp = UniDecAPIMethods.InputMethods.SetupInputs();


			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();
			unsafe
			{
				fixed (float* xarrayPtr = &xarray[0], yarrayPtr = &yarray[0])
				{
					inp.dataMZ = xarrayPtr;
					inp.dataInt = yarrayPtr;
					inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);

					int result = AlgorithmTesting.TestSetStartEnds(inp, config);
					Assert.AreEqual(10994, result);
				}
			}
		}
		[Test]
		public void TestMemoryObjectAllocationToHeapConfigPtr()
		{
			// initialize the config and inp structs for data transfer to and from C
			Config config = new();
			InputUnsafe inp = new();

			// sets the config to deafult
			config = UniDecAPIMethods.ConfigMethods.ModifyConfigToDefault(config);
			// modifies some MsDataScan-dependent field of config
			UniDecAPIMethods.ConfigMethods.SetupConfig(ref config, scan);
			// initializes the inp object to deafult. 			
			inp = UniDecAPIMethods.InputMethods.SetupInputs();

			// gets the m/z and int arrays from msdatascan. 
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();

			unsafe
			{   // initializes pointers to the arrays. 
				fixed (float* xarrayPtr = &xarray[0], yarrayPtr = &yarray[0])
				{
					// assigns the pointers to the inp.dataMZ and inp.dataInt fields in the inp struct
					inp.dataMZ = xarrayPtr;
					inp.dataInt = yarrayPtr;
					// function calculates some values that are data dependent
					inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);

					// actual test
					int result = AlgorithmTesting.MemoryObjectAllocationToHeapConfigPtr(config, inp);
					Assert.AreEqual(1, result);
				}
			}
		}
		[Test]
		public void TestFreeDecon()
		{
			int result = AlgorithmTesting.TestFreeDecon();
			Assert.AreEqual(1, result);
		}
		[Test]
		public void TestSetupAndAllocateMemoryToDecon()
		{
			int result = AlgorithmTesting.TestSetupAndAllocateMemoryToDecon();
			Assert.AreEqual(1, result);
		}
		[Test]
		public void TestSetupAndReturnDecon()
		{
			// error was caused by returning a pointer to a managed type (in this case float[])
			Decon decon = AlgorithmTesting.TestSetupAndReturnDecon();
			PrintProperties(decon);
		}
		[Test]
		public void TestUniDecDeconvolution()
		{

		}
		[Test]
		public unsafe void TestBarr()
		{
			// setup the config struct
			UniDecAPIMethods.ConfigMethods.CreateAndSetupConfig(scan, out Config config);

			// setup the input struct
			InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();

			// assign inp the x and y array data
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();

			int numberElementsInBarr = config.lengthmz * config.numz;
			char[] barr = new char[numberElementsInBarr];

			fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
			{
				inp.dataInt = ptrYarray;
				inp.dataMZ = ptrXarray;
				inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
				byte[] barrResult = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.barr, numberElementsInBarr);
				Assert.AreEqual(0, barr[0]);
				Assert.AreEqual(numberElementsInBarr, barrResult.Length); 
			}
		}
		[Test]
		public unsafe void TestReadInputsByValue()
		{
			// setup the config struct
			UniDecAPIMethods.ConfigMethods.CreateAndSetupConfig(scan, out Config config);

			// setup the input struct
			InputUnsafe inp = UniDecAPIMethods.InputMethods.SetupInputs();

			// assign inp the x and y array data
			float[] xarray = scan.MassSpectrum.XArray.ConvertDoubleArrayToFloat();
			float[] yarray = scan.MassSpectrum.YArray.ConvertDoubleArrayToFloat();

			int numberElementsInBarr = config.lengthmz * config.numz;
			char[] barr = new char[numberElementsInBarr];

			fixed (float* ptrXarray = &xarray[0], ptrYarray = &yarray[0])
			{
				inp.dataInt = ptrYarray;
				inp.dataMZ = ptrXarray;
				inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
			}
		}

		[Test]
		public unsafe void TestCreateAndSetupConfig()
		{
			UniDecAPIMethods.ConfigMethods.CreateAndSetupConfig(scan, out Config config);
			PrintProperties(config);
		}

		[Test]
		public unsafe void TestPtrToArrayConversions()
		{
			float[] testFloatArray = new float[] { 1F, 2F, 3F };
			char[] testCharArray = new char[] { '0', '1', '0' }; 
			
			fixed(float* ptrFloatArray = &testFloatArray[0])
			{
				float[] resultFloatArray = UniDecAPIMethods.UtilityMethods.PtrToArray(ptrFloatArray, testFloatArray.Length);
				Assert.AreEqual(testFloatArray, resultFloatArray); 
			}
			
			fixed(char* ptrCharArray = &testCharArray[0])
			{
				char[] resultCharArray = UniDecAPIMethods.UtilityMethods.PtrToArray(ptrCharArray, testCharArray.Length);
				Assert.AreEqual(testCharArray, resultCharArray); 
			}
		}
		[Test]
		public unsafe void TestCharArrayMarshallingToC()
		{

			char[] testCharArray = {'0', '0', '0'};
			char[] expectedArray = { '1', '1', '1' };
			byte[] testByteArray = { 0, 0, 0 }; 
			TestingCharArrayMarshalling(testCharArray, 3);
			TestingCharArrayMarshalling(testByteArray, 3); 
			// this test fails because char is 2 byte in C#, but one byte in C. 
			Assert.AreNotEqual(expectedArray, testCharArray);
			// This test passses because byte is one byte in C# and char is one byte in C. 
			Assert.AreEqual(expectedArray, testByteArray); 
		}

		[DllImport("TestDLL.dll", EntryPoint = "TestingCharArrayMarshalling")]
		private static unsafe extern void _TestingCharArrayMarshalling(char* array, int arrayLength);
		[DllImport("TestDLL.dll", EntryPoint = "TestingCharArrayMarshalling")]
		private static unsafe extern void _TestingCharArrayMarshalling(byte* array, int arrayLength);
		public static unsafe void TestingCharArrayMarshalling(char[] array, int arrayLength)
		{
			fixed (char* charPtr = &array[0])
			{
				_TestingCharArrayMarshalling(charPtr, arrayLength);
			}
		}
		public static unsafe void TestingCharArrayMarshalling(byte[] array, int arrayLength)
		{
			fixed (byte* bytePtr = &array[0])
			{
				_TestingCharArrayMarshalling(bytePtr, arrayLength);
			}
		}
		[Test]
		public void TestSetupAndMakeIsotopes()
		{


		}


	}
}
