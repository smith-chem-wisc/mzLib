using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static unsafe class Convolution
	{
		public static int SetStartEnds(Config config, ref InputUnsafe inp, ref int[] starttab, ref int[] endtab, float threshold)
		{
			// changed InputUnsafe to pass by reference here, but I think I could've also passed it directly and 
			// created the pointer with a fixed statement as well. Probably better to pass by reference though. 

			int maxlength = 1;
			for (int i = 0; i < config.lengthmz; i++)
			{
				float point = inp.dataMZ[i] - threshold;
				int start, end;
				if (point < inp.dataMZ[0] && config.speedyflag == 0)
				{
					//start = (int)((point - inp->dataMZ[0]) / (inp->dataMZ[1] - inp->dataMZ[0]));
					start = 0 - Nearfast(inp.dataMZ, (float)2 * inp.dataMZ[0] - point, config.lengthmz);
				}
				else
				{
					start = Nearfast(inp.dataMZ, point, config.lengthmz);
				}
				starttab[i] = start;

				point = inp.dataMZ[i] + threshold;
				if (point > inp.dataMZ[config.lengthmz - 1] && config.speedyflag == 0)
				{
					//end = config.lengthmz - 1 + (int)((point - inp->dataMZ[config.lengthmz - 1]) / (inp->dataMZ[config.lengthmz - 1] - inp->dataMZ[config.lengthmz - 2]));
					end = config.lengthmz - 1 + Nearfast(inp.dataMZ, (float)2 * inp.dataMZ[0] - point, config.lengthmz);
				}
				else
				{
					end = Nearfast(inp.dataMZ, point, config.lengthmz);
				}
				endtab[i] = end;
				if (end - start > maxlength) { maxlength = end - start; }
				//printf("%d %d\n", start, end);
			}
			//printf("Max Length: %d\t", maxlength);
			return maxlength;
		}
		public static int Nearfast(float* dataMz, float point, int numdat)
		{
			return _Nearfast(dataMz, point, numdat);
		}
		public static void DeconvolveBaseline(Config config, InputUnsafe inp, Decon decon)
		{
			fixed (float* baselinePtr = &decon.baseline[0])
			{
				_DeconvolveBaseline(config.lengthmz, inp.dataMZ, inp.dataInt, baselinePtr, Math.Abs(config.mzsig));
			}

		}
		public static float Reconvolve(int lengthmz, int numz, int maxlength, int[] starttab, int[] endtab,
			float[] mzdist, float[] blur, float[] newblur, int speedyflag, byte[] barr)
		{
			fixed (int* starttabPtr = &starttab[0], endtabPtr = &endtab[0])
			{
				fixed (float* mzdistPtr = &mzdist[0], blurPtr = &blur[0], newblurPtr = &newblur[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						return _Reconvolve(lengthmz, numz, maxlength, starttabPtr, endtabPtr,
							mzdistPtr, blurPtr, newblurPtr, speedyflag, barrPtr);
					}
				}
			}

		}

		// Private methods 
		[DllImport("TestDLL.dll", EntryPoint = "nearfast")]
		private static extern int _Nearfast(float* dataMz, float point, int numdat);

		[DllImport("TestDLL.dll", EntryPoint = "deconvolve_baseline")]
		private static extern void _DeconvolveBaseline(int lengthmz, float* dataMZ, float* dataInt, float* baseline,
			float mzsig);

		[DllImport("TestDLL.dll", EntryPoint = "Reconvolve")]
		private static extern float _Reconvolve(int lengthmz, int numz, int maxlength, int* starttab,
			int* endtab, float* mzdist, float* blur, float* newblur, int speedyflag, byte* barr);

	}
}
