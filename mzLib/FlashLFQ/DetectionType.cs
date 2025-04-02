namespace FlashLFQ
{

    public enum DetectionType
    {
        MSMS,                           // The peak is detected from MS2ID
        MBR,                            // The peak is detected from MBR
        MSMSAmbiguousPeakfinding,       // The peak is detected from more than one MS2ID
        IsoTrack_MBR,                   // The peak is detected from MBR by IsoTrack
        IsoTrack_Ambiguous,             // Multiple peptides are mapped to this peak by IsoTracker.
        MSMSIdentifiedButNotQuantified, // We have MS2ID but no peak for quantification
        NotDetected,                    // We don't have MS2ID either peak for quantification
    }
}