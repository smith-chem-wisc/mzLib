namespace MzLibSpectralAveraging
{
    /// <summary>
    ///  Enum of options for spectrum binning. Likely that SpectrumBinning will be deprecated.
    /// Use MrsNoiseEstimate as the default. 
    /// </summary>
    public enum SpectrumMergingType
    {
        SpectrumBinning,
        MostSimilarSpectrum, 
        MrsNoiseEstimate
    }
}
