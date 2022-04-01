using System;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static unsafe class Isotopes
	{
		public static float Isotopemid(float mass, float* isoparams)
		{
			float a, b, c;
			a = isoparams[4];
			b = isoparams[5];
			c = isoparams[6];
			return a + b * (float)Math.Pow(mass, c);
		}

		public static float Isotopesig(float mass, float* isoparams)
		{
			float a, b, c;
			a = isoparams[7];
			b = isoparams[8];
			c = isoparams[9];
			return a + b * (float)Math.Pow(mass, c);
		}

		public static float Isotopealpha(float mass, float* isoparams)
		{
			float a, b;
			a = isoparams[0];
			b = isoparams[1];
			return a * (float)Math.Exp(-mass * b);
		}

		public static float Isotopebeta(float mass, float* isoparams)
		{
			float a, b;
			a = isoparams[2];
			b = isoparams[3];
			return a * (float)Math.Exp(-mass * b);
		}
		[DllImport("TestDLL.dll", EntryPoint = "SetupAndMakeIsotopes")]
		private static extern void _SetupAndMakeIsotopes(Config config, InputUnsafe inp);
		public static void SetupAndMakeIsotopes(Config config, InputUnsafe inp)
		{
			_SetupAndMakeIsotopes(config, inp);
		}
		[DllImport("TestDLL.dll", EntryPoint = "setup_isotopes")]
		private static extern IsotopeStruct _SetupIsotopes(float* isoparams, int* isotopepos, float* isotopeval,
			float* mtab, int* ztab, byte* barr, float* dataMZ, int lengthmz, int numz);
		public static IsotopeStruct SetupIsotopes(Config config, InputUnsafe inp)
		{
			fixed (float* isoparamsPtr = &InputUnsafe.isoparams[0])
			{
				return _SetupIsotopes(isoparamsPtr, inp.isotopeops, inp.isotopeval, inp.mtab, inp.nztab,
				inp.barr, inp.dataMZ, config.lengthmz, config.numz);
			}
		}
		[DllImport("TestDLL.dll", EntryPoint = "make_isotopes")]
		private static extern void _MakeIsotopes(float* isoparams, int* isotopepos, float* isotopeval, float* mtab, int* ztab,
			byte* barr, float* dataMZ, int lengthmz, int numz, float minmid, float maxmid, float maxsig);
		public static void MakeIsotopesFromAPI(Config config, InputUnsafe inp, IsotopeStruct isoStruct)
		{
			fixed (float* isoparamsPtr = &InputUnsafe.isoparams[0])
			{
				_MakeIsotopes(isoparamsPtr, inp.isotopeops, inp.isotopeval, inp.mtab, inp.nztab, inp.barr, inp.dataMZ,
				config.lengthmz, config.numz, isoStruct.minmid, isoStruct.maxmid, isoStruct.maxsig);
			}
		}
		public static void MakeIsotopes(Config config, InputUnsafe inp, IsotopeStruct isos)
		{
			fixed (float* isoparamsPtr = &InputUnsafe.isoparams[0])
			{
				float massdiff = 1.0026f;
				int isostart = 0;
				int isoend = (int)(isos.maxmid + 4 * isos.maxsig);

				int isolength = isoend - isostart;

				float[] isorange = new float[isolength];
				int[] isoindex = new int[isolength];

				for (int i = 0; i < isolength; i++)
				{
					isorange[i] = (isostart + i) * massdiff;
					isoindex[i] = (isostart + i);
				}

				for (int i = 0; i < config.lengthmz; i++)
				{
					for (int j = 0; j < config.numz; j++)
					{
						if (inp.barr[UniDecAPIMethods.UtilityMethods.index2D(config.numz, i, j)] - 48 == 1)
						{
							float mz = inp.dataMZ[i];
							int z = inp.nztab[i];

							// need to declare variables outside fixed if going to use outside of fixed. 
							float alpha, amp, beta, tot, mid, sig;
							float mass = inp.mtab[ArrayIndexing.Index2D(config.numz, i, j)];


							mid = Isotopemid(mass, isoparamsPtr);
							sig = Isotopesig(mass, isoparamsPtr);

							alpha = Isotopealpha(mass, isoparamsPtr);
							amp = (1.0f - alpha) / (sig * 2.50662827f);
							beta = Isotopebeta(mass, isoparamsPtr);
							tot = 0;


							for (int k = 0; k < isolength; k++)
							{
								float newmz = mz + (isorange[k] + (float)z);
								int pos = Convolution.Nearfast(inp.dataMZ, newmz, config.lengthmz);

								float e = alpha * (float)Math.Exp(-isoindex[k] * beta);
								float g = amp * (float)Math.Exp(-Math.Pow(isoindex[k] - mid, 2) / (2 * Math.Pow(sig, 2)));
								float temp = e + g;
								tot += temp;
								if (tot > 0)
								{
									temp *= 1 / tot;
								}
								inp.isotopeval[ArrayIndexing.Index3D(config.numz, isolength, i, j, k)] = temp;
							}
						}
					}
				}
			}
		}

		public static void MonoisotopicToAverageMass(Config config, InputUnsafe inp, Decon decon, byte[] barr)
		{
			float[] newblur = new float[config.lengthmz * config.numz];
			for (int i = 0; i < config.lengthmz; i++)
			{
				for (int j = 0; j < config.numz; j++)
				{
					if (barr[ArrayIndexing.Index2D(config.numz, i, j)] == 1)
					{
						float topval = decon.blur[ArrayIndexing.Index2D(config.numz, i, j)];
						for (int k = 0; k < config.isolength; k++)
						{
							int pos = inp.isotopeops[ArrayIndexing.Index3D(config.numz, config.isolength, i, j, k)];
							float val = inp.isotopeval[ArrayIndexing.Index3D(config.numz, config.isolength, i, j, k)];
							newblur[ArrayIndexing.Index2D(config.numz, pos, j)] += topval * val;
						}
					}
				}
			}
		}
	}
}
