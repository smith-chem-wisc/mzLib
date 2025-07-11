namespace MassSpectrometry
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
        /// Refers to the index of the scan in an array that contains scans that were indexed 
        /// (i.e., for indexing MS1 peaks, the array only contains MS1 scans)
        /// </summary>
        public int ZeroBasedScanIndex { get; }
        /// <summary>
        /// Represents the mass, either as m/z or a neutral mass, depending on the context
        /// </summary>
        public double M { get; }
    }
}
