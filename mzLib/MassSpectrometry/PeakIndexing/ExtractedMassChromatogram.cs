using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.PeakIndexing
{
    public class ExtractedMassChromatogram : ExtractedIonChromatogram
    {
        public ExtractedMassChromatogram(List<IIndexedPeak> peaks) : base(peaks)
        {
            
        }

    }
}
