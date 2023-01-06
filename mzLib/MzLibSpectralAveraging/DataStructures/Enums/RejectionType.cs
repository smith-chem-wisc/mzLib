namespace MzLibSpectralAveraging
{
    /// <summary>
    /// List of rejection types implemented. 
    /// </summary>
    public enum RejectionType
    {
        NoRejection,
        MinMaxClipping,
        PercentileClipping,
        SigmaClipping,
        WinsorizedSigmaClipping,
        AveragedSigmaClipping,
        BelowThresholdRejection,
        Thermo
    }


}
