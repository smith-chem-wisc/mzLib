using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MassSpectrometry
{
    public class ExtractedIonChromatogram
    {
        private List<Datum> XicData;

        public ExtractedIonChromatogram(List<Datum> xicData)
        {
            this.XicData = xicData;
        }

        /// <summary>
        /// Returns the XIC data, ordered by retention time. 
        /// Each datum's X field is the m/z value, Y is the intensity, and Z is the retention time.
        /// </summary>
        public IEnumerable<Datum> Data
        {
            // return the data, ordered by retention time
            get { return XicData.OrderBy(p => p.X).AsEnumerable(); }
        }
    }
}
