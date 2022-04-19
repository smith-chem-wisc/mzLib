using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public static unsafe class FitFunctions
	{
		public static void KillB(float* I, byte[] B, float intthresh, int lengthmz, int numz,
			int isolength, int* isotopepos, float* isotopeval)
		{
			fixed (byte* bPtr = &B[0])
			{
				_KillB(I, bPtr, intthresh, lengthmz, numz, isolength,
					isotopepos, isotopeval);
			}
		}
		public static void KillB(float* I, byte* B, float intthresh, int lengthmz, int numz,
			int isolength, int* isotopepos, float* isotopeval)
		{
			_KillB(I, B, intthresh, lengthmz, numz, isolength, isotopepos, isotopeval);
		}
		public static void KillB(InputUnsafe inp, Config config, byte[] barr)
		{
			fixed (byte* barrPtr = &barr[0])
			{
				_KillB(inp.dataInt, barrPtr, config.intthresh, config.lengthmz, config.numz,
					config.isolength, inp.isotopeops, inp.isotopeval);
			}
		}
		[DllImport("TestDLL.dll", EntryPoint = "KillB")]
		private static extern void _KillB(float* I, byte* b, float intthresh, int lengthmz, int numz,
			int isolength, int* isotopepos, float* isotopeval);
	}
}
