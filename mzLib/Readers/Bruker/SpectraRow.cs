using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers
{
    // class used as a one-to-one map of sqlite object in bruker analysis.baf file for conversion to MsDataScan. 
    internal class SpectraTableRow
	{
		public int Id { get; set; }
		public double Rt { get; set; }
		public int Segment { get; set; }
		public int AcquisitionKey { get; set; }
		public int? Parent { get; set; }
		public int MzAcqRangeLower { get; set; }
		public int MzAcqRangeUpper { get; set; }
		public double SumIntensity { get; set; }
		public double MaxIntensity { get; set; }
		public long? TransformatorId { get; set; }
		public long? ProfileMzId { get; set; }

		public long? ProfileIntensityId { get; set; }
		public long? LineIndexId { get; set; }
		public long? LineMzId { get; set; }
		public long? LineIntensityId { get; set; }
		public long? LineIndexWidthId { get; set; }
		public long? LinePeakAreaId { get; set; }
		public long? LineSnrId { get; set; }
	}
}
