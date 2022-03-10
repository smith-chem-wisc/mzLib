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
		[MarshalAs(UnmanagedType.LPArray)]
		public float[] dataMZ;
		[MarshalAs(UnmanagedType.LPArray)]
		public float[] dataInt;
		[MarshalAs(UnmanagedType.LPArray)]
		public float[] testmasses;
		[MarshalAs(UnmanagedType.LPArray)]
		public int[] nztab;
		[MarshalAs(UnmanagedType.LPArray)]
		public float[] mtab;
		[MarshalAs(UnmanagedType.LPTStr)] // check to see if this is compatible
		public string barr;
		[MarshalAs(UnmanagedType.LPArray)]
		public int[] isotopeops;
		[MarshalAs(UnmanagedType.LPArray)]
		public float[] isotopeval;
		[MarshalAs(UnmanagedType.LPArray)]
		public float[] isoparams; 
	}
}
