using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ.Interfaces
{
    /// <summary>
    /// Represent an object that can be traced over a separation domain
    /// E.g., A chromatographic peak trace composed of multiple IsotopicEnvelopes
    /// </summary>
    public interface ITraceable<T> where T : ISingleScanDatum
    {
        /// <summary>
        /// The most intense point in the trace
        /// </summary>
        public T Apex { get; }
        /// <summary>
        /// A list of data points that compose the ITraceable object
        /// This list must be ordered by the separation domain in ascending order! (e.g., retention time)
        /// </summary>
        public List<T> ScanOrderedPoints { get; }

        public void CalculateIntensityForThisFeature(bool integrate);

        public void SetSplitLocation(T splitPoint);
    }
}
