using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    /// <summary>
    /// An ISingleScanDatum represents information that exists in a single mass spectrometric scan! 
    /// E.g., a single m/z peak, a single isotopic envelope
    /// </summary>
    public interface ISingleScanDatum
    {
        public double Mz { get; }
        public double Intensity { get; }

        /// <summary>
        /// The position of the datum in the separation domain,
        /// e.g., retention time, ion mobility
        /// </summary>
        public double RelativeSeparationValue { get; }
        public int ZeroBasedScanIndex { get; }
    }
}
