using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static unsafe class MathUtilities
	{
		[DllImport("TestDLL.dll", EntryPoint = "Max")]
		private static extern float _Max(float* blur, int length);
		public static float Max(float[] blur, int length)
		{
			fixed (float* blurPtr = &blur[0])
			{
				return _Max(blurPtr, length);
			}
		}
		public static float Max(float* blur, int length)
		{
			return _Max(blur, length);
		}
		public static void AplyCutoff1D(float* array, float cutoff, int lengthmz)
		{
			for (int i = 0; i < lengthmz; i++)
			{
				if (array[i] < cutoff)
				{
					array[i] = 0;
				}
			}
		}
	}
}
