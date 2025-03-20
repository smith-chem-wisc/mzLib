namespace FlashLFQ
{
    public enum DetectionType
    {
        MSMS,
        MBR,
        IsoTrack_MBR, // MBR detected by IsoTrack
        IsoTrack_Ambiguous, // Ambiguous(more than two Id in one peak) detected by IsoTrack
        NotDetected,
        MSMSAmbiguousPeakfinding,
        MSMSIdentifiedButNotQuantified,
        Imputed,
        Default // Default value, will be removed in the future
    }
}