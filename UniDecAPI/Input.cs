using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	[StructLayout(LayoutKind.Sequential)]
	public struct Input
	{
		public float[] dataMZ;
		public float[] dataInt;
		public float[] testmasses;
		public int[] nztab;
		public float[] mtab;
		public string barr;
		public int[] isotopeops;
		public float[] isotopeval;
		public float[] isoparams; 
	}
}
