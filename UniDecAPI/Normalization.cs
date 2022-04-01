using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public unsafe static class Normalization
	{
		public static void SimpNormSum(int length, float[] data)
		{
			fixed (float* dataPtr = &data[0])
			{
				simp_norm_sum(length, dataPtr);
			}
		}

		[DllImport("TestDLL.dll", EntryPoint = "simp_norm_sum")]
		private static extern void simp_norm_sum(int length, float* data);

	}
}
