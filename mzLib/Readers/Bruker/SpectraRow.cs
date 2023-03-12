using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Readers.Bruker
{
	internal class SpectraTableRow
	{
		public int Id { get; set; }
		public double Rt { get; set; }
		public int Segment { get; set; }
		public int AcquisitionKey { get; set; }
		public int? ParentSpectrum { get; set; }
		public int MzAcqRangeLower { get; set; }
		public int MzAcqRangeUpper { get; set; }
		public double SumIntensity { get; set; }
		public double MaxIntensity { get; set; }
		public int? TransformatorId { get; set; }
		public int? ProfileMzId { get; set; }

		public int? ProfileIntensityId { get; set; }
		public int? LineIndexId { get; set; }
		public int? LineMzId { get; set; }
		public int? LineIntensityId { get; set; }
		public int? LineIndexWidthId { get; set; }
		public int? LinePeakAreaId { get; set; }
		public int? LineSnrId { get; set; }
	}
}
