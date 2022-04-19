using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	[StructLayout(LayoutKind.Sequential)]
	public struct Decon
	{
		public float[] fitdat;
		public float[] baseline;
		public float[] noise;
		public float[] massgrid;
		public float[] massaxis;
		public float[] massaxisval;
		public float[] blur;
		public float[] newblur;
		public float[] peakx;
		public float[] peaky;
		public float[] dscores;
		public float error;
		public float rsquared;
		public float iterations;
		public float uniscore;
		public float conv;
		public float threshold;
		public int plen;
		public int mlen; 
	}
	public unsafe struct DeconUnsafe
	{
		public float* fitdat;
		public float* baseline;
		public float* noise;
		public float* massgrid;
		public float* massaxis;
		public float* massaxisval;
		public float* blur;
		public float* newblur;
		public float* peakx;
		public float* peaky;
		public float* dscores;
		public float error;
		public float rsquared;
		public float iterations;
		public float uniscore;
		public float conv;
		public float threshold;
		public int plen;
		public int mlen;
	}
	public class DeconResults
	{
		public float[] Fitdat { get; private set; }
		public float[] Baseline { get; private set; }
		public float[] Noise { get; private set; }
		public float[] MassGrid { get; private set; }
		public float[] MassAxis { get; private set; }
		public float[] MassAxisVal { get; private set; }
		public float[] PeakX { get; private set; }
		public float[] PeakY { get; private set; }
		public float[] DScores { get; private set; }
		private float Iterations { get; set; }
		public float UniScore { get; private set; }

		public DeconResults(DeconUnsafe decon, Config config)
		{
			unsafe
			{
				Fitdat = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.fitdat, config.lengthmz);
				Baseline = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.baseline, config.lengthmz * config.numz);
				Noise = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.noise, config.lengthmz * config.numz);
				MassGrid = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.massgrid, decon.mlen * config.numz);
				MassAxis = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.massaxis, decon.mlen);
				MassAxisVal = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.massaxisval, decon.mlen);
				PeakX = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.peakx, decon.mlen);
				PeakY = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.peaky, decon.mlen);
				DScores = UniDecAPIMethods.UtilityMethods.PtrToArray(decon.dscores, decon.plen);
				UniScore = decon.uniscore; 
			}
		}
	}
}
