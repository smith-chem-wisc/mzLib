using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices; 

namespace UniDecAPI
{
	[StructLayout(LayoutKind.Sequential)]
	public struct IsotopeStruct
	{
		public float minmid;
		public float minsig;
		public float maxmid;
		public float maxsig;
		int isolength; 
	}
}
