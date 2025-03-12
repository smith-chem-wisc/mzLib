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
    public interface IIndexedPeak
    {
        public double Intensity { get; }

        /// <summary>
        /// The position of the datum in the separation domain,
        /// e.g., retention time, ion mobility
        /// </summary>
        public double RetentionTime { get; }
        public int ZeroBasedMs1ScanIndex { get; }
        public double Mz { get; }
    }
}
