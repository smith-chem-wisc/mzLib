using System.Runtime.InteropServices;
using System; 

namespace UniDecAPI
{
	public static unsafe class Scoring
	{
		[DllImport("TestDLL.dll", EntryPoint = "performScoring")]
		// UniScore method in C requires that decon is passed as a pointer and the function returns a modified version of 
		// the decon struct. So we need the "In" keyword to specify that a pointer needs to get passed and the out 
		// to make sure we're getting the modified Decon struct out. 
		private static extern float _UniScore(Config config, [In,Out] Decon decon, InputUnsafe inp, float scoreThreshold);
		public static float UniScore(Config config, Decon decon, InputUnsafe inp, float scoreThreshold)
		{
			return _UniScore(config, decon, inp, scoreThreshold);
		}
		public static float UniScorePorted(Config config, Decon decon, InputUnsafe inp, float threshold, int windowSize)
		{
			decon.peakx = new float[decon.mlen];
			decon.peaky = new float[decon.mlen];
			int plen = 0;

			fixed (float* peakxPtr = &decon.peakx[0], peakyPtr = &decon.peaky[0])
			{
				plen = PeakDetect(inp.dataMZ, inp.dataInt, config.lengthmz, windowSize,
					threshold, peakxPtr, peakyPtr);
			}
			decon.plen = plen;
			
			// once plen is calculated, reinitialize peakx and peaky to a float[plen]. 
			decon.peakx = new float[plen];
			decon.peaky = new float[plen];
			decon.dscores = new float[plen];

			PeakNorm(decon, config);
			return ScoreFromPeaks(decon, config, inp, threshold);
		}

		public static int PeakDetect(Decon decon, Config config)
		{
				return _PeakDetect(decon.massaxis, decon.massaxisval, decon.mlen,
							config.peakwin, config.peakthresh, decon.peakx, decon.peaky);

		}
		public static void PeakNorm(Decon decon, Config config)
		{
			_PeakNorm(decon.peaky, decon.plen, config.peaknorm); 
		}
		public static float ScoreFromPeaks(Decon decon, Config config, InputUnsafe inp, float threshold)
		{
			return _ScoreFromPeaks(decon.plen, decon.peaky, decon.dscores, config, decon,
				inp, threshold); 
		}
		public unsafe static int IsPeak(float* dataMz, float* dataInt, int lengthmz, float window, 
			float thresh, int index)
		{
			float xval = dataMz[index];
			float yval = dataInt[index]; 

			if(yval < thresh) { return 0; }
			for(int i = 0; i < lengthmz; i++)
			{
				float temp = dataMz[i]; 
				if(Math.Abs(temp - xval) <= window)
				{
					float tempy = dataInt[i]; 
					if(tempy > yval)
					{
						return 0; 
					}
					if (tempy == yval & i < index)
					{
						return 0; 
					}
				}
			}
			return 1; 
		}
		public unsafe static int PeakDetect(float* dataMz, float* dataInt, int lengthmz, float window, 
			float thresh, float* peakx, float* peaky)
		{
			int plen = 0;
			float max = MathUtilities.Max(dataInt, lengthmz); 
			for(int i = 0; i < lengthmz; i++)
			{
				if (IsPeak(dataMz, dataInt, lengthmz, window, thresh * max, i) == 1)
				{
					peakx[plen] = dataMz[i];
					peaky[plen] = dataInt[i];
					plen++; 
				}
			}
			return plen; 
		}

		
		public static unsafe float UScore(Config config, float* dataMz, float* dataInt, float* mzgrid,
			int* nztab, float mlow, float mhigh, float peak)
		{
			return _Uscore(config, dataMz, dataInt, mzgrid, nztab, mlow, mhigh, peak); 
		}
		public static void GetFWHMS(Config config, int plen, int mlen, float* massaxis, float* massaxisval, 
			float* peakx, out float[] fwhmlow, out float[] fwhmhigh, out float[] badfwhm)
		{
			fwhmlow = new float[plen];
			fwhmhigh = new float[plen];
			badfwhm = new float[plen];

			fixed(float* fwhmlowPtr = &fwhmlow[0], fwhmhighPtr = &fwhmhigh[0], badfwhmPtr = &badfwhm[0])
			{
				_GetFWHMS(config, plen, mlen, massaxis, massaxisval, peakx, fwhmlowPtr, fwhmhighPtr, badfwhmPtr);
				fwhmlow = UniDecAPIMethods.UtilityMethods.PtrToArray(fwhmlowPtr, plen);
				fwhmhigh = UniDecAPIMethods.UtilityMethods.PtrToArray(fwhmhighPtr, plen);
				badfwhm = UniDecAPIMethods.UtilityMethods.PtrToArray(badfwhmPtr, plen);
			}

		}
		public static float MScore(Config config, int mlen, float* massaxis, float* masssum, 
			float* massgrid, float mlow, float mhigh, float peak)
		{
			return _Mscore(config, mlen, massaxis, masssum, massgrid, mlow, mhigh, peak); 
		}
		public static float CsScore(Config config, int mlen, float* massaxis, float* masssum, 
			float* massgrid, float mlow, float mhigh, float peak)
		{
			return _Csscore(config, mlen, massaxis, masssum, massgrid, mlow, mhigh, peak); 
		}
		public static float FScore(Config config, int plen, int mlen, float* massaxis, float* masssum, 
			float* peakx, float height, float mlow, float mhigh, float peak, int badfwhm)
		{
			return _Fscore(config, plen, mlen, massaxis, masssum, peakx, height, mlow, mhigh, peak, badfwhm); 
		}

		public static float ScoreFromPeaksPorted(int plen, float* peakx, float* peaky,
			float* dscores, Config config, Decon decon, InputUnsafe inp, float threshold)
		{
			float xfwhm = 2;
			float[] fwhmHigh;
			float[] fwhmLow;
			float[] badFwhm;

			fixed (float* massaxisPtr = &decon.massaxis[0], massaxisvalPtr = &decon.massaxisval[0],
				peakxPtr = &decon.peakx[0], newblurPtr = &decon.newblur[0], massgridPtr = &decon.massgrid[0])
			{
				GetFWHMS(config, plen, decon.mlen, massaxisPtr, massaxisvalPtr, peakxPtr, 
					out fwhmLow, out fwhmHigh, out badFwhm);

				float numerator = 0f;
				float denominator = 0f;
				float uniscore = 0f;

				for (int i = 0; i < plen; i++)
				{

					float m = peakx[i];
					float ival = peaky[i];
					float l = m - (m - fwhmLow[i]) * xfwhm;
					float h = m + (fwhmHigh[i] - m) * xfwhm;
					int index = Convolution.Nearfast(massaxisPtr, m, decon.mlen);
					float height = massaxisvalPtr[index];

					float usc = UScore(config, inp.dataMZ, inp.dataInt, newblurPtr, inp.nztab, l, h, m);
					float msc = MScore(config, decon.mlen, massaxisPtr, massaxisvalPtr, massgridPtr, l, h, m); 
					float cssc = CsScore(config, decon.mlen, massaxisPtr, massaxisvalPtr, massgridPtr, l,h,m);
					float fsc = FScore(config, plen, decon.mlen, massaxisPtr, massaxisvalPtr, peakxPtr, height, 
						fwhmLow[i], fwhmHigh[i], m, (int)badFwhm[i]); // not sure why badFwhm isn't int*. Explicit cast now before I fix it later. 

					float dsc = usc * msc * cssc * fsc;
					dscores[i] = dsc; 
					if(dsc > threshold)
					{
						numerator += ival * ival * dsc;
						denominator += ival * ival; 
					}
				}
				if(denominator != 0)
				{
					uniscore = decon.rsquared * numerator / denominator; 
				}
				return uniscore; 
			}
		}

		[DllImport("TestDLL.dll", EntryPoint = "peak_detect")]
		private static extern int _PeakDetect(float* massaxis, float* massaxisval, int mlen,
			float peakwin, float peakthresh, float* peakx, float* peaky);

		[DllImport("TestDLL.dll", EntryPoint = "peak_detect")]
		private static extern int _PeakDetect([In, Out] float[] massaxis, [In, Out] float[] massaxisval, int mlen,
			float peakwin, float peakthresh, [In, Out] float[] peakx, [In, Out] float[] peaky);

		[DllImport("TestDLL.dll", EntryPoint = "peak_norm")]
		private static extern void _PeakNorm([In, Out] float[] peaky, int plen, int peaknorm);

		[DllImport("TestDLL.dll", EntryPoint = "score_from_peaks")]
		private static extern float _ScoreFromPeaks(int plen, [In, Out] float[] peaky, [In, Out] float[] dscores,
			Config config, [In, Out] Decon decon, InputUnsafe inp, float threshold);

		[DllImport("TestDLL.dll", EntryPoint = "uscore")]
		private static extern float _Uscore(Config config, float* dataMz, float* dataInt, float* mzgrid,
			int* nztab, float mlow, float mhigh, float peak);

		[DllImport("TestDLL.dll", EntryPoint = "mscore")]
		private static extern float _Mscore(Config config, int mlen, float* massaxis, float* masssum,
			float* massgrid, float mlow, float mhigh, float peak);

		[DllImport("TestDLL.dll", EntryPoint = "cssscore")]
		private static extern float _Csscore(Config config, int mlen, float* massaxis, float* masssum,
			float* massgrid, float mlow, float mhigh, float peak);

		[DllImport("TestDLL.dll", EntryPoint = "fscore")]
		private static extern float _Fscore(Config config, int plen, int mlen, float* massaxis,
			float* masssum, float* peakx, float height, float mlow, float mhigh, float peak, int badfwhm);

		[DllImport("TestDll.dll", EntryPoint = "get_fwhms")]
		private static extern void _GetFWHMS(Config config, int plen, int mlen, float* massaxis,
			float* massaxisval, float* peakx, float* fwhmlow, float* fwhmhigh, float* badfwhm);


	}


}
