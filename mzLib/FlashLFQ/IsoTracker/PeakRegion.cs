
using System;
namespace FlashLFQ.IsoTracker
{
    /// <summary>  
    /// An PeakRegion should contain Apex Retention Time, Start and End Retention Time in the region.  
    /// </summary>  
    public class PeakRegion : IEquatable<PeakRegion>
    {
        public double ApexRT { get; private set; }
        public double StartRT { get; private set; }
        public double EndRT { get; private set; }
        /// <summary>  
        /// This property calculates the width of the peak region by subtracting the start retention time from the end retention time.  
        /// </summary>  
        public double Width { get { return EndRT - StartRT; } }

        public PeakRegion(double apexRT, double startRT, double endRT)
        {
            this.ApexRT = apexRT;
            this.StartRT = startRT;
            this.EndRT = endRT;
        }

        public bool Equals(PeakRegion other)
        {
            return ApexRT == other.ApexRT && StartRT == other.StartRT && EndRT == other.EndRT;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(ApexRT, StartRT, EndRT);
        }
    }
}
