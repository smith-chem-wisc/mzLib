using System;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public partial class UniDecAPIMethods
	{
		public static class DeconMethods
		{
			[DllImport("UniDecMinimal.dll", EntryPoint = "SetupDecon")]
			private static extern Decon _SetUpDecon(); 
			public static Decon SetupDecon()
			{
				return _SetUpDecon(); 
			}
			[DllImport("UniDecMinimal.dll", EntryPoint = "FreeDecon")]
			private static extern void _FreeDecon(Decon decon); 
			public static void FreeDecon(Decon decon)
			{
				_FreeDecon(decon); 
			}
			[DllImport("UniDecMinimal.dll", EntryPoint = "Average")]
			private static extern float _Average(int length, float[] xarray); 
			public static float Average(int length, float[] xarray)
			{
				return _Average(length, xarray); 
			}
		}
	}
}
