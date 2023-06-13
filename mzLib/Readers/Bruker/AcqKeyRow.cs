using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
	// class used as a one-to-one map of sqlite object in bruker analysis.baf file for conversion to MsDataScan. 
	internal class AcqKeyRow
	{
		public int Id { get; set; }
		// 0 = positive, 1 = negative 
		public int Polarity { get; set; }
		// 0 = MS
		// 2 = MS/MS
		// 4 = in-source CID
		// 5 = broadband CID
		// 255 = unknown
		public int ScanMode { get; set; }
		// 1 = (axial or orthogonal) TOF, linear detection mode
		// 2 = (axial or orthogonal) TOF, reflector detection mode
		// 33 = FTMS
		// 255 = unknown
		public int AcquisitionMode { get; set; }
		// 0 = MS (no fragmentation),
		// 1 = MS/MS
		public int MsLevel { get; set; }
	}
}
