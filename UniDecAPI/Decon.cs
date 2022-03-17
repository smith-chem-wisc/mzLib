using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	public unsafe struct Decon
	{
		public float* fitdat;
		public float* baseline;
		public float* noise;
		public float* massgride;
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
}
