using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using System.Runtime.InteropServices;

namespace UniDecAPI
{
	public static class UniDecDeconvolution
	{ 
		public static void DeconvoluteWithUniDec(this MzSpectrum spectrum, out DeconResults deconResults)
		{
			UniDecDeconEngine deconEngine = new();
			deconResults = deconEngine.Run(spectrum); 
		}
	}
	public class UniDecDeconEngine 
	{
		Config config;
		DeconUnsafe decon;
		InputUnsafe inp;
		UnmanagedHandling unHandler;
		byte[] barr; 

		public DeconResults Run(MzSpectrum spectrum)
		{
			DeconResults results; 
			using (unHandler = new UnmanagedHandling())
			{
				Setup(spectrum);
				Processing();
				results = ScoreAndOutput(); 
			}
			return results; 
		}

		private void Setup(MzSpectrum spectrum)
		{
			unsafe
			{
				// create the config, input and decon structs. 
				UniDecAPIMethods.ConfigMethods.CreateAndSetupConfig(spectrum, out config);
				inp = UniDecAPIMethods.InputMethods.SetupInputs();
				decon = new DeconUnsafe();

				// copy the x- and y-axis to a pointer and assign it to the input struct. 
				inp.dataInt = (float*)unHandler.AllocateToUnmanaged(spectrum.YArray.ConvertDoubleArrayToFloat());
				inp.dataMZ = (float*)unHandler.AllocateToUnmanaged(spectrum.XArray.ConvertDoubleArrayToFloat());
				
				// further inp setup after the data from the spectrum is filled. 
				inp = UniDecAPIMethods.InputMethods.ReadInputsByValue(inp, config);
				UniDecAPIMethods.UtilityMethods.SetLimits(config, inp);
				IsotopeStruct isoStruct = Isotopes.SetupIsotopes(config, inp);
				config.isolength = Math.Abs(isoStruct.isolength);

				inp.isotopeops = (int*)unHandler.AllocateToUnmanaged(config.isolength * config.lengthmz * config.numz, typeof(int));
				inp.isotopeval = (float*)unHandler.AllocateToUnmanaged(config.isolength * config.lengthmz * config.numz, typeof(float));

				Isotopes.MakeIsotopesFromAPI(config, inp, isoStruct);

				// I believe barr stands for "Bad array". 0 if bad; 1 if good. 
				int numberOfElementsInBarr = config.lengthmz * config.numz;

				barr = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.barr, numberOfElementsInBarr);
				// C is returning a char, so this function needs to convert the char to the equivalent C# byte. 
				UniDecAPIMethods.UtilityMethods.ConvertASCIIBytesFromCToByteInCS(ref barr);
			}
		}
		private unsafe void Processing()
		{
			float threshold = config.psthresh * Math.Abs(config.mzsig) * config.peakshapeinflate;

			int* starttab = (int*)unHandler.AllocateToUnmanaged(config.lengthmz, typeof(int));
			int* endtab = (int*)unHandler.AllocateToUnmanaged(config.lengthmz, typeof(int));

			int maxlength = Convolution.SetStartEnds(config, ref inp, starttab, endtab, threshold);

			int pslen = config.lengthmz * maxlength;
			float* mzdist = (float*)unHandler.AllocateToUnmanaged(pslen, typeof(float));
			float* rmzdist = (float*)unHandler.AllocateToUnmanaged(pslen, typeof(int));
			int makereverse = 1;
			// makepeakshape2d is very slow (>20s) 
			MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse, starttab, endtab, maxlength);

			int zlength = 1 + 2 * (int)config.zsig;
			int mlength = 1 + 2 * (int)config.msig;
			int* mind = (int*)unHandler.AllocateToUnmanaged(mlength, typeof(int));
			float* mdist = (float*)unHandler.AllocateToUnmanaged(mlength, typeof(int));
			for (int i = 0; i < mlength; i++)
			{
				mind[i] = i - (mlength - 1) / 2;
				if (config.msig != 0)
				{
					mdist[i] = (float)(Math.Exp(-Math.Pow((i - (zlength - 1) / 2.0), 2)) / (2.0 * config.zsig * config.zsig));
				}
				else
				{
					mdist[i] = 1;
				}
			}
			int* zind = stackalloc int[zlength];
			float* zdist = stackalloc float[zlength];

			int numclose = mlength * zlength;
			int* closemind = (int*)unHandler.AllocateToUnmanaged(numclose, typeof(int));
			int* closezind = (int*)unHandler.AllocateToUnmanaged(numclose, typeof(int));
			float* closeval = (float*)unHandler.AllocateToUnmanaged(numclose, typeof(int));
			int* closeind = (int*)unHandler.AllocateToUnmanaged(numclose * config.lengthmz * config.numz, typeof(int));
			float* closearray = (float*)unHandler.AllocateToUnmanaged(numclose * config.lengthmz * config.numz, typeof(int));

			for (int k = 0; k < numclose; k++)
			{
				closemind[k] = mind[k % mlength];
				closezind[k] = zind[(int)k / mlength];
				closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength];
			}

			Normalization.simp_norm_sum(mlength, mdist);
			Normalization.simp_norm_sum(zlength, zdist);
			Normalization.simp_norm_sum(numclose, closeval);

			Blur._MakeSparseBlur(numclose, inp.barr, closezind, closemind, inp.mtab, inp.nztab,
				inp.dataMZ, closeind, closeval, closearray, config);

			int badness = 1;
			for (int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if (barr[i] == 1)
				{
					badness = 0;
				}
			}
			if (badness == 1)
			{
				throw new InvalidOperationException("Badness = 1...");
			}

			float dmax = MathUtilities.Max(inp.dataInt, config.lengthmz);
			float betafactor = 1;
			if (dmax > 1)
			{
				betafactor = dmax;
			}

			UniDecAPIMethods.UtilityMethods.KillB_CS(inp.dataInt, ref barr, config.intthresh, config.lengthmz, config.numz, config.isolength,
				inp.isotopeops, inp.isotopeval);
			// DirectUniDecPort.FitFunctions.KillB(inp.dataInt, barr, 
			//	config.intthresh, config.lengthmz, config.numz, 
			//	config.isolength, inp.isotopeops, inp.isotopeval); 

			decon.blur = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			decon.newblur = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			float* oldblur = stackalloc float[config.lengthmz * config.numz];
			decon.baseline = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));
			decon.noise = (float*)unHandler.AllocateToUnmanaged(config.lengthmz * config.numz, typeof(float));

			for (int i = 0; i < config.lengthmz; i++)
			{
				float val = inp.dataInt[i] / ((float)(config.numz + 2));
				if (config.baselineflag == 1)
				{

					decon.baseline[i] = val;
					decon.noise[i] = val;
				}

				for (int j = 0; j < config.numz; j++)
				{
					if (barr[ArrayIndexing.Index2D(config.numz, i, j)] == 1)
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

			// copy decon.blur to oldblur: 
			oldblur = decon.blur;
			// copy decon.newblur to decon.blur: 
			decon.newblur = decon.blur;

			float[] dataInt2 = UniDecAPIMethods.UtilityMethods.PtrToArray(inp.dataInt, config.lengthmz);


			Convolution.DeconvolveBaseline(config, inp, decon);
			float blurmax = 0F;
			decon.conv = 0F;

			fixed (float* dataInt2Ptr = &dataInt2[0])
			{
				Blur.PerformIterations(ref decon, config, inp, betafactor, maxlength,
					starttab, endtab, mzdist, numclose, closeind, closearray, zlength, mlength,
					closemind, closezind, mdist, dataInt2Ptr, zdist, barr, rmzdist, oldblur);
			}


			if (config.peakshapeinflate != 1 && config.mzsig != 0)
			{
				if (config.speedyflag == 0)
				{
					MZPeak.MakePeakShape2D(config, inp, mzdist, rmzdist, makereverse: 0, starttab, endtab, maxlength);
				}
				else
				{
					MZPeak.MakePeakShape1D(config, inp, threshold, mzdist, rmzdist, makeReverse: 0);
				}
			}

			blurmax = MathUtilities.Max(decon.blur, config.lengthmz * config.numz);
			float cutoff = 0F;
			if (blurmax != 0)
			{
				cutoff = 0.000001F;
			}

			ArrayIndexing.ApplyCutoff1D(decon.blur, blurmax * cutoff, config.lengthmz * config.numz);

			decon.fitdat = (float*)unHandler.AllocateToUnmanaged(config.lengthmz, typeof(float));


			decon.error = ErrorFunctions.ErrFunctionSpeedy(config, ref decon, barr, inp.dataInt, maxlength, inp.isotopeops, inp.isotopeval,
				starttab, endtab, mzdist);
			//decon.error = ErrorFunctions.ErrFunctionSpeedy(config, decon, barr, inp.dataInt, maxlength,
			//	inp.isotopeops, inp.isotopeval, starttab, endtab, mzdist);

			if (config.intthresh != -1)
			{
				for (int i = 0; i < config.lengthmz - 1; i++)
				{
					if (inp.dataInt[i] == 0 && inp.dataInt[i + 1] == 0)
					{
						decon.fitdat[i] = 0F;
						decon.fitdat[i + 1] = 0F;
					}
				}
			}
			// not tested yet. 
			if (config.isotopemode == 2)
			{
				Isotopes.MonoisotopicToAverageMass(config, inp, decon, barr);
			}

			float newblurmax = blurmax;
			if (config.rawflag == 0 || config.rawflag == 2)
			{
				if (config.mzsig != 0)
				{
					newblurmax = Convolution._Reconvolve(config.lengthmz, config.numz, maxlength,
						starttab, endtab, mzdist, decon.blur, decon.newblur, config.speedyflag, inp.barr);
				}
				else
				{
					decon.newblur = decon.blur;
				}
			}
			float massmax = config.masslb;
			float massmin = config.massub;

			if (config.fixedmassaxis == 0)
			{
				for (int i = 0; i < config.lengthmz; i++)
				{
					for (int j = 0; j < config.numz; j++)
					{
						if (decon.newblur[ArrayIndexing.Index2D(config.numz, i, j)] * barr[ArrayIndexing.Index2D(config.numz, i, j)] > newblurmax * cutoff)
						{
							float testmax = inp.mtab[ArrayIndexing.Index2D(config.numz, i, j)] + threshold * inp.nztab[j] + config.massbins;
							float testmin = inp.mtab[ArrayIndexing.Index2D(config.numz, i, j)] - threshold * inp.nztab[j];

							//To prevent really wierd decimals
							testmin = (float)Math.Round(testmin / config.massbins) * config.massbins;
							testmax = (float)Math.Round(testmax / config.massbins) * config.massbins;

							if (testmax > massmax) { massmax = testmax; }
							if (testmin < massmin) { massmin = testmin; }
						}
					}
				}
			}
			else { massmax = config.massub; massmin = config.masslb; }

			//Checks to make sure the mass axis is good and makes a dummy axis if not
			decon.mlen = (int)((massmax - massmin) / config.massbins);
			if (decon.mlen < 1)
			{
				massmax = config.massub;
				massmin = config.masslb;
				decon.mlen = (int)((massmax - massmin) / config.massbins);

				//Declare the memory

				decon.massaxis = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massaxisval = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massgrid = (float*)unHandler.AllocateToUnmanaged(decon.mlen * config.numz, typeof(float));

				//Create the mass axis
				for (int i = 0; i < decon.mlen; i++)
				{
					decon.massaxis[i] = massmin + i * config.massbins;
				}
				decon.uniscore = 0;
			}
			else
			{

				//Declare the memory
				decon.massaxis = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massaxisval = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
				decon.massgrid = (float*)unHandler.AllocateToUnmanaged(decon.mlen * config.numz, typeof(float));


				//Create the mass axis
				for (int i = 0; i < decon.mlen; i++)
				{
					decon.massaxis[i] = massmin + i * config.massbins;
				}
			}

			MassIntensityDetermination.IntegrateMassIntensities(config, ref decon, inp);
		}
		private unsafe DeconResults ScoreAndOutput()
		{
			float scorethresh = 0f;
			config.peakwin = 20;

			decon.peakx = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));
			decon.peaky = (float*)unHandler.AllocateToUnmanaged(decon.mlen, typeof(float));

			decon.uniscore = Scoring.UniScorePorted(config, ref decon, inp, scorethresh, config.peakwin, unHandler);

			return new DeconResults(decon, config); 
		}

	}

}
