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
		public fixed float isoparams[10];
	}
}
