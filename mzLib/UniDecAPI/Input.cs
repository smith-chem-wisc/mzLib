using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	[StructLayout(LayoutKind.Sequential)]
	public unsafe struct InputUnsafe
	{
		public float* dataMZ;
		public float* dataInt;
		public float* testmasses;
		public int* nztab;
		public float* mtab;
		public byte* barr;
		public int* isotopeops;
		public float* isotopeval;
		public static float[] isoparams = {
			1.00840852e+00F,
			1.25318718e-03F,
			2.37226341e+00F,
			8.19178000e-04F,
			-4.37741951e-01F,
			6.64992972e-04F,
			9.94230511e-01F,
			4.64975237e-01F,
			1.00529041e-02F,
			5.81240792e-01F
		};
	}
}
