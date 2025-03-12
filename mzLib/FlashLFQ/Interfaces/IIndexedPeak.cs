using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    /// <summary>
    /// An IIndexedPeak represents information that exists in a single mass spectrometric scan 
    /// E.g., a single m/z peak
    /// </summary>
    public interface IIndexedPeak
    {
        public double Intensity { get; }
        public double RetentionTime { get; }
        /// <summary>
        /// Explicitly refers to the index of the MS1 scan in an array that contains ONLY MS1 scans, ordered by retention time
        /// </summary>
        public int ZeroBasedMs1ScanIndex { get; }
        public double Mz { get; }
    }
}
