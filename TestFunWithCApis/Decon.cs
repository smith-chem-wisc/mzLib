using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestUniDec
{
	public partial class UniDecClasses
	{
		public struct Decon
		{
			float[] fitdat;
			float[] baseline;
			float[] noise;
			float[] massgrid;
			float[] massaxis;
			float[] massaxisval;
			float[] blur;
			float[] newblur;
			float[] peakx;
			float[] peaky;
			float[] dscores;
			float error;
			float rsquared;
			int iterations;
			float uniscore;
			float conv;
			float threshold;
			int mlen;
			int plen; 
		}

	}
}
