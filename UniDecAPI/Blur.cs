using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static unsafe class Blur
	{
		public static void MakeSparseBlur(int numclose, byte* barr, int* closezind,
			int* closemind, float* mtab, int* nztab, float* dataMZ, int* closeind, float* closeval,
			float* closearray, Config config)
		{
			_MakeSparseBlur(numclose, barr, closezind, closemind, mtab,
				nztab, dataMZ, closeind, closeval, closearray, config);
		}
		public static void MakeSparseBlur(InputUnsafe inp, Config config, int numclose, byte[] barr,
			int[] closezind, int[] closemind, int[] closeind,
			float[] closeval, float[] closearray)
		{
			fixed (byte* barrPtr = &barr[0])
			{
				fixed (int* closezindPtr = &closezind[0], closemindPtr = &closemind[0],
					closeindPtr = &closeind[0])
				{
					fixed (float* closevalPtr = &closeval[0], closearrayPtr = &closearray[0])
					{
						_MakeSparseBlur(numclose, barrPtr, closezindPtr, closemindPtr, inp.mtab, inp.nztab,
							inp.dataMZ, closeindPtr, closevalPtr, closearrayPtr, config);
					}
				}
			}
		}
		public static void CreateInitialBlur(Decon decon, InputUnsafe inp, Config config)
		{
			for (int i = 0; i < config.lengthmz; i++)
			{
				float val = inp.dataInt[i] / (float)(config.numz + 2);
				if (config.baselineflag == 1)
				{

					decon.baseline[i] = val;
					decon.noise[i] = val;
				}

				for (int j = 0; j < config.numz; j++)
				{
					if (inp.barr[ArrayIndexing.Index2D(config.numz, i, j)] == 1)
					{
						if (config.isotopemode == 0)
						{
							decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = val;
						}
						else { decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 1; }
					}
					else
					{
						decon.blur[ArrayIndexing.Index2D(config.numz, i, j)] = 0;
					}
				}
			}
		}
		public static void SoftargmaxTransposed(Decon decon, Config config, InputUnsafe inp,
			float betaFactor, byte[] barr, int maxlength, int[] starttab, int[] endtab, float[] rmzdist)
		{
			float beta = Math.Abs(config.beta / betaFactor);
			fixed (float* rmzdistPtr = &rmzdist[0], blurPtr = &decon.blur[0])
			{
				fixed (int* starttabPtr = &starttab[0], endtabPtr = &endtab[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						_SoftArgmaxTransposed(blurPtr, config.lengthmz, config.numz, beta,
							barrPtr, maxlength, config.isolength, inp.isotopeops, inp.isotopeval, config.speedyflag,
							starttabPtr, endtabPtr, rmzdistPtr, config.mzsig);
					}
				}
			}

		}
		public static void Softargmax(float[] blur, int lengthmz, int numz, float beta)
		{
			fixed (float* blurPtr = &blur[0])
			{
				_Softargmax(blurPtr, lengthmz, numz, beta);
			}

		}
		public static void PointSmoothing(float[] blur, byte[] barr, int lengthmz, int numz, int width)
		{
			fixed (float* blurPtr = &blur[0])
			{
				fixed (byte* barrPtr = &barr[0])
				{
					_PointSmoothing(blurPtr, barrPtr, lengthmz, numz, width);
				}
			}

		}
		public static void PointSmoothingPeakWidth(int lengthmz, int numz, int maxlength,
			int[] starttab, int[] endtab, float[] mzdist, float[] blur, int speedyflag, byte[] barr)
		{
			fixed (float* mzdistPtr = &mzdist[0], blurPtr = &blur[0])
			{
				fixed (int* starttabPtr = &starttab[0], endtabPtr = &endtab[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						_PointSmoothingPeakWidth(lengthmz, numz, maxlength, starttabPtr, endtabPtr,
							mzdistPtr, blurPtr, speedyflag, barrPtr);
					}
				}
			}
		}
		public static void BlurItMean(Config config, Decon decon, int numclose, int[] closeind,
			byte[] barr, float[] closeArray)
		{
			fixed (float* newblurPtr = &decon.newblur[0], closearrayPtr = &closeArray[0],
				blurPtr = &decon.blur[0])
			{
				fixed (int* closeIndPtr = &closeind[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						_BlurItMean(config.lengthmz, config.numz, numclose, closeIndPtr,
							newblurPtr, blurPtr, barrPtr, closearrayPtr, config.zerolog);
					}
				}
			}

		}
		public static void BlurItHybrid1(int lengthmz, int numz, int zlength, int mlength,
			int[] closeind, int[] closemind, int[] closezind, float[] mdist, float[] zdist, float[] newblur,
			float[] blur, byte[] barr, float[] closearray, float zerolog)
		{
			fixed (int* closeindPtr = &closeind[0], closemindPtr = &closemind[0], closezindPtr = &closezind[0])
			{
				fixed (float* mdistPtr = &mdist[0], zdistPtr = &zdist[0], newblurPtr = &newblur[0],
					blurPtr = &blur[0], closearrayPtr = &closearray[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						_BlurItHybrid1(lengthmz, numz, zlength, mlength, closeindPtr, closemindPtr,
							closezindPtr, mdistPtr, zdistPtr, newblurPtr, blurPtr, barrPtr, closearrayPtr, zerolog);
					}
				}
			}

		}
		public static void BlurIt(int lengthmz, int numz, int numclose, int[] closeind,
			float[] closearray, float[] newblur, float[] blur, byte[] barr)
		{
			fixed (int* closeindPtr = &closeind[0])
			{
				fixed (float* closearrayPtr = &closearray[0], newblurPtr = &newblur[0], blurPtr = &blur[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						_BlurIt(lengthmz, numz, numclose, closeindPtr, closearrayPtr, newblurPtr, blurPtr, barrPtr);
					}
				}
			}

		}
		public static void PerformIterations(ref Decon decon, Config config, InputUnsafe inp, float betafactor, int maxlength,
			int[] starttab, int[] endtab, float[] mzdist, int numclose, int[] closeind, float[] closearray, int zlength, int mlength,
			int[] closemind, int[] closezind, float[] mdist, float[] dataInt2, float[] zdist, byte[] barr, float[] rmzdist, float[] oldblur)
		{
			int off = 0;
			for (int iterations = 0; iterations < Math.Abs(config.numit); iterations++)
			{
				decon.iterations = iterations;
				if (config.beta > 0 && iterations > 0)
				{

					Softargmax(decon.blur, config.lengthmz, config.numz, config.beta / betafactor);
					//printf("Beta %f\n", beta);
				}
				else if (config.beta < 0 && iterations > 0)
				{
					SoftargmaxTransposed(decon, config, inp, betafactor, barr, maxlength, starttab, endtab, rmzdist);
				}

				if (config.psig >= 1 && iterations > 0)
				{
					PointSmoothing(decon.blur, barr, config.lengthmz, config.numz, Math.Abs((int)config.psig));
					//printf("Point Smoothed %f\n", config.psig);
				}
				else if (config.psig < 0 && iterations > 0)
				{
					PointSmoothingPeakWidth(config.lengthmz, config.numz, maxlength, starttab, endtab, mzdist, decon.blur, config.speedyflag, barr);
				}


				//Run Blurs
				if (config.zsig >= 0 && config.msig >= 0)
				{
					BlurItMean(config, decon, numclose, closeind, barr, closearray);
				}
				else if (config.zsig > 0 && config.msig < 0)
				{
					BlurItHybrid1(config.lengthmz, config.numz, zlength, mlength, closeind, closemind,
						closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
				}
				else if (config.zsig < 0 && config.msig > 0)
				{
					BlurItHybrid2(config.lengthmz, config.numz, zlength, mlength,
						closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
				}
				else
				{
					BlurIt(config.lengthmz, config.numz, numclose, closeind, closearray, decon.newblur, decon.blur, barr);
				}

				//Run Richardson-Lucy Deconvolution
				DeconvolveIterationSpeedy(config.lengthmz, config.numz, maxlength,
					decon.newblur, decon.blur, barr, config.aggressiveflag, dataInt2,
					config.isolength, inp.isotopeops, inp.isotopeval, starttab, endtab, mzdist, rmzdist, config.speedyflag,
					config.baselineflag, decon.baseline, decon.noise, config.mzsig, inp.dataMZ, config.filterwidth, config.psig);

				//Determine the metrics for conversion. Only do this every 10% to speed up.
				if ((config.numit < 10 || iterations % 10 == 0 || iterations % 10 == 1 || iterations > 0.9 * config.numit))
				{
					float diff = 0;
					float tot = 0;
					for (int i = 0; i < config.lengthmz * config.numz; i++)
					{
						if (barr[i] == 1)
						{
							diff += (float)Math.Pow(((double)decon.blur[i] - (double)oldblur[i]), 2);
							tot += decon.blur[i];
						}
					}
					if (tot != 0) { decon.conv = (diff / tot); }
					else { decon.conv = 12345678F; }

					//printf("Iteration: %d Convergence: %f\n", iterations, decon.conv);
					if (decon.conv < 0.000001F)
					{
						if (off == 1 && config.numit > 0)
						{

							break;
						}
						off = 1;
					}
				}
				oldblur = decon.blur;
			}

		}
		[DllImport("TestDLL.dll", EntryPoint = "MakeSparseBlur")]
		private static extern void _MakeSparseBlur(int numclose, byte* barr, int* closezind,
			int* closemind, float* mtab, int* nztab, float* dataMZ, int* closeind, float* closeval,
			float* closearray, Config config);

		[DllImport("TestDLL.dll", EntryPoint = "softargmax")]
		private static extern void _Softargmax(float* blur, int lengthmz, int numz, float beta);

		[DllImport("TestDLL.dll", EntryPoint = "softargmax_transposed")]
		private static extern void _SoftArgmaxTransposed(float* blur, int lengthmz, int numz, float beta,
			byte* barr, int maxlength, int isolength, int* isotopepos, float* isotopeval, int speedyflag,
			int* starttab, int* endtab, float* mzdist, float mzsig);

		[DllImport("TestDLL.dll", EntryPoint = "point_smoothing")]
		private static extern void _PointSmoothing(float* blur, byte* barr, int lengthmz, int numz, int width);

		[DllImport("TestDLL.dll", EntryPoint = "point_smoothing_peak_width")]
		private static extern void _PointSmoothingPeakWidth(int lengthmz, int numz, int maxlength, int* starttab, int* endtab, float* mzdist, float* blur, int speedyflag, byte* barr);

		[DllImport("TestDLL.dll", EntryPoint = "blur_it_mean")]
		private static extern void _BlurItMean(int lengthmz, int numz, int numclose, int* closeind, float* newblur,
			float* blur, byte* barr, float* closearray, float zerolog);

		[DllImport("TestDLL.dll", EntryPoint = "blur_it_hybrid1")]
		private static extern void _BlurItHybrid1(int lengthmz, int numz, int zlength, int mlength,
			int* closeind, int* closemind, int* closezind, float* mdist, float* zdist, float* newblur,
			float* blur, byte* barr, float* closearray, float zerolog);

		[DllImport("TestDLL.dll", EntryPoint = "blur_it_hybrid2")]
		private static extern void _BlurItHybrid2(int lengthmz, int numz, int zlength, int mlength,
			int* closeind, int* closemind, int* closezind, float* mdist, float* zdist, float* newblur,
			float* blur, byte* barr, float* closearray, float zerolog);
		public static void BlurItHybrid2(int lengthmz, int numz, int zlength, int mlength,
			int[] closeind, int[] closemind, int[] closezind, float[] mdist, float[] zdist, float[] newblur,
			float[] blur, byte[] barr, float[] closearray, float zerolog)
		{
			fixed (int* closeindPtr = &closeind[0], closemindPtr = &closemind[0], closezindPtr = &closezind[0])
			{
				fixed (float* mdistPtr = &mdist[0], zdistPtr = &zdist[0], newblurPtr = &newblur[0],
					blurPtr = &blur[0], closearrayPtr = &closearray[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						_BlurItHybrid2(lengthmz, numz, zlength, mlength, closeindPtr, closemindPtr, closezindPtr,
							mdistPtr, zdistPtr, newblurPtr, blurPtr, barrPtr, closearrayPtr, zerolog);
					}
				}
			}
		}
		[DllImport("TestDLL.dll", EntryPoint = "blur_it")]
		private static extern void _BlurIt(int lengthmz, int numz, int numclose, int* closeind,
			float* closearray, float* newblur, float* blur, byte* barr);

		[DllImport("TestDLL.dll", EntryPoint = "deconvolve_iteration_speedy")]
		private static extern void _DeconvolveIterationSpeedy(int lengthmz, int numz, int maxlength,
			float* blur, float* blur2, byte* barr, int aggressiveflag, float* dataInt,
			int isolength, int* isotopepos, float* isotopeval, int* starttab, int* endtab,
			float* mzdist, float* rmzdist, int speedyflag, int baselineflag, float* baseline,
			float* noise, float mzsig, float* dataMZ, float filterwidth, float psig);

		public static void DeconvolveIterationSpeedy(int lengthmz, int numz, int maxlength,
			float[] blur, float[] blur2, byte[] barr, int aggressiveflag, float[] dataInt,
			int isolength, int* isotopepos, float* isotopeval, int[] starttab, int[] endtab,
			float[] mzdist, float[] rmzdist, int speedyflag, int baselineflag, float[] baseline,
			float[] noise, float mzsig, float* dataMZ, float filterwidth, float psig)
		{
			fixed (float* blurPtr = &blur[0], blur2Ptr = &blur2[0], dataIntPtr = &dataInt[0],
				mzdistPtr = &mzdist[0], rmzdistPtr = &rmzdist[0], baselinePtr = &baseline[0], noisePtr = &noise[0])
			{
				fixed (int* starttabPtr = &starttab[0], endtabPtr = &endtab[0])
				{
					fixed (byte* barrPtr = &barr[0])
					{
						_DeconvolveIterationSpeedy(lengthmz, numz, maxlength,
							blurPtr, blur2Ptr, barrPtr, aggressiveflag, dataIntPtr, isolength, isotopepos,
							isotopeval, starttabPtr, endtabPtr, mzdistPtr, rmzdistPtr, speedyflag,
							baselineflag, baselinePtr, noisePtr, mzsig, dataMZ, filterwidth, psig);
					}
				}
			}

		}


	}

}
