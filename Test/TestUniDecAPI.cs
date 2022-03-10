using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using UniDecAPI; 

namespace Test
{
	class TestUniDecAPI
	{
		[Test]
		public void TestDefaultConfig()
		{
			Config testConfig = UniDecAPIMethods.ConfigMethods.SetDefaultConfig();
			PrintProperties(testConfig); 
		}
		[Test]
		public void TestDecon()
		{
			Decon decon = new(); 
			try
			{
				decon = UniDecAPIMethods.DeconMethods.SetupDecon();
				PrintProperties(decon); 
			}
			finally
			{
				UniDecAPIMethods.DeconMethods.FreeDecon(decon); 
			}
		}
		[Test]
		public void TestInput()
		{

		}
		[Test]
		public void TestAverage()
		{
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
		
	}
}
