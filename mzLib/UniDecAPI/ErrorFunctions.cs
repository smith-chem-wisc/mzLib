using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static unsafe class ErrorFunctions
	{
		[DllImport("TestDLL.dll", EntryPoint = "errfunspeedy")]
		private static extern float _ErrFunctionSpeedy(Config config, Decon decon, byte* barr, float* dataInt, int maxlength,
			int* isotopepos, float* isotopeval, int* starttab, int* endtab, float* mzdist, ref float rsquared,
			float* fitdat, float* blur, float* baseline);
		[DllImport("TestDLL.dll", EntryPoint = "errfunspeedy")]
		private static extern float _ErrFunctionSpeedy(Config config, ref DeconUnsafe decon, byte* barr, float* dataInt, int maxlength,
			int* isotopepos, float* isotopeval, int* starttab, int* endtab, float* mzdist, ref float rsquared,
			float* fitdat, float* blur, float* baseline);

		public static float ErrFunctionSpeedy(Config config, Decon decon, byte[] barr, float* dataInt, int maxlength,
			int* isotopepos, float* isotopeval, int[] starttab, int[] endtab, float[] mzdist)
		{
			// the issue is that it requires quite a few calls to Decon's arrays. So I need to actually break it out 
			// into many more pointers. 
			fixed (int* starttabPtr = &starttab[0], endtabPtr = &endtab[0])
			{
				fixed (float* mzdistPtr = &mzdist[0], fitdatPtr = &decon.fitdat[0],
					blurPtr = &decon.blur[0], baselinePtr = &decon.baseline[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						return _ErrFunctionSpeedy(config, decon, barrPtr, dataInt, maxlength, isotopepos, isotopeval,
									starttabPtr, endtabPtr, mzdistPtr, ref decon.rsquared, fitdatPtr, blurPtr, baselinePtr);
					}
				}
			}
		}
		public static float ErrFunctionSpeedy(Config config, ref DeconUnsafe decon, byte[] barr, float* dataInt, int maxlength,
			int* isotopepos, float* isotopeval, int* starttab, int* endtab, float* mzdist)
		{
			fixed (byte* barrPtr = &barr[0])
			{
				return _ErrFunctionSpeedy(config, ref decon, barrPtr, dataInt, maxlength, isotopepos, isotopeval,
							starttab, endtab, mzdist, ref decon.rsquared, decon.fitdat, decon.blur, decon.baseline);
			}
		}
	}
}
