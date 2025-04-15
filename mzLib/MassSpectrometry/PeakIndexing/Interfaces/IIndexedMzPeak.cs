namespace MassSpectrometry
{
    /// <summary>
    /// An IIndexedMzPeak represents information that exists in a single mass spectrometric scan 
    /// E.g., a single m/z peak
    /// </summary>
    public interface IIndexedMzPeak
    {
        public double Intensity { get; }
        public double RetentionTime { get; }
        /// <summary>
        /// Refers to the index of the scan in an array that contains scans that were indexed 
        /// (i.e., for indexing MS1 peaks, the array only contains MS1 scans)
        /// </summary>
        public int ZeroBasedScanIndex { get; }
        public double Mz { get; }
    }
}
