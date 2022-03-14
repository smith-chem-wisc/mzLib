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

namespace Test
{
	// minimize your unsafe code calling!!! 
	// If you have to use it, minimize the unsafe scope!! 

	class TestUniDecAPI
	{
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
			foreach(var field in o.GetType().GetFields())
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
		[TestCase("LowResMS1ScanForDecon.raw")]
		public void TestPerformUniDecDeconvolution(string file)
		{
			var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
			List<MsDataScan> testScan = ThermoRawFileReader.LoadAllStaticData(path).GetAllScansList();
			var scan = testScan[0];

			// Decon deconResult = scan.PerformUniDecDeconvolution(0, 0); 
			scan.TestPerformUniDecDecon(); 
		}
	}
}
