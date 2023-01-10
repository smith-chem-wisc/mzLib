namespace MzLibSpectralAveraging
{
    /// <summary>
    /// List of rejection types implemented. 
    /// </summary>
    public enum OutlierRejectionType
    {
        NoRejection,
        MinMaxClipping,
        PercentileClipping,
        SigmaClipping,
        WinsorizedSigmaClipping,
        AveragedSigmaClipping,
        BelowThresholdRejection,
    }


}
