using System;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public partial class UniDecAPIMethods
	{
		public static class DeconMethods
		{
			[DllImport("TestDLL.dll", EntryPoint = "SetupDecon")]
			private static extern Decon _SetUpDecon(); 
			public static Decon SetupDecon()
			{
				return _SetUpDecon(); 
			}
			[DllImport("TestDLL.dll", EntryPoint = "FreeDecon")]
			private static extern void _FreeDecon(Decon decon); 
			public static void FreeDecon(Decon decon)
			{
				_FreeDecon(decon); 
			}
			[DllImport("TestDLL.dll", EntryPoint = "Average")]
			private static extern float _Average(int length, float[] xarray); 
			public static float Average(int length, float[] xarray)
			{
				return _Average(length, xarray); 
			}

			[DllImport("TestDLL.dll", EntryPoint = "MainDeconvolution")]
			private static extern Decon _MainDeconvolution(Config config, InputUnsafe inp, int silent, int verbose);
			public static Decon MainDeconvolution(Config config, InputUnsafe inp, int silent, int verbose)
			{
				return _MainDeconvolution(config, inp, silent, verbose); 
			}
			[DllImport("TestDLL.dll", EntryPoint = "TestingMainDeconvolution")]
			private static extern int _TestingMainDeconvolution();
			public static int TestingMainDeconvolution()
			{
				return _TestingMainDeconvolution(); 
			}
		}
	}
}
