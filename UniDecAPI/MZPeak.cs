using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static unsafe class MZPeak
	{
		public static void MakePeakShape2D(Config config, InputUnsafe inp, float[] mzdist, float[] rmzdist,
			 int makereverse, int[] starttab, int[] endtab, int maxlength)
		{
			fixed (float* mzdistPtr = &mzdist[0], rmzdistPtr = &rmzdist[0])
			{
				fixed (int* starttabPtr = &starttab[0], endtabPtr = &endtab[0])
				{
					float mzsig = Math.Abs(config.mzsig) * config.peakshapeinflate;
					_MakePeakShape2D(config.lengthmz, maxlength, starttabPtr, endtabPtr,
						inp.dataMZ, mzsig, config.psfun, config.speedyflag, mzdistPtr, rmzdistPtr, makereverse);
				}
			}
		}
		public static void MakePeakShape2D(Config config, InputUnsafe inp, float* mzdist, float* rmzdist,
			 int makereverse, int* starttab, int* endtab, int maxlength)
		{
					float mzsig = Math.Abs(config.mzsig) * config.peakshapeinflate;
					_MakePeakShape2D(config.lengthmz, maxlength, starttab, endtab,
						inp.dataMZ, mzsig, config.psfun, config.speedyflag, mzdist, rmzdist, makereverse);
		}
		public static void MakePeakShape1D(Config config, InputUnsafe inp, float threshold, float[] mzdist,
			float[] rmzdist, int makeReverse)
		{
			fixed (float* mzdistPtr = &mzdist[0], rmzdistPtr = &rmzdist[0])
			{
				_MakePeakShape1D(inp.dataMZ, threshold, config.lengthmz, config.speedyflag,
					Math.Abs(config.mzsig), config.psfun, mzdistPtr, rmzdistPtr, makeReverse);
			}

		}
		public static void MakePeakShape1D(Config config, InputUnsafe inp, float threshold, float* mzdist,
			float* rmzdist, int makeReverse)
		{
			_MakePeakShape1D(inp.dataMZ, threshold, config.lengthmz, config.speedyflag,
				Math.Abs(config.mzsig), config.psfun, mzdist, rmzdist, makeReverse);
		}

		[DllImport("TestDLL.dll", EntryPoint = "MakePeakShape2D")]
		private static extern void _MakePeakShape2D(int lengthmz, int maxlength,
			int* starttab, int* endtab, float* dataMZ, float mzsig, int psfun, int speedyflag,
			float* mzdist, float* rmzdist, int makereverse);

		[DllImport("TestDLL.dll", EntryPoint = "MakePeakShape1D")]
		private static extern void _MakePeakShape1D(float* dataMZ, float threshold, int lengthmz, int speedyflag, float mzsig,
			int psfun, float* mzdist, float* rmzdist, int makereverse);
	}
}
