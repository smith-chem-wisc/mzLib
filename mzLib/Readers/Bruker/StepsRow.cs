using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IO.BrukerFileReader
{
	internal class StepsRow
	{
		public int TargetSpectrum { get; set; }
		public int Number { get; set; }
		public int IsolationType { get; set; }
		public int ReactionType { get; set; }
		// 0 = MS; 1 = MS/MS
		public int MsLevel { get; set; }
		public double ParentMass { get; set; }
	}
}
