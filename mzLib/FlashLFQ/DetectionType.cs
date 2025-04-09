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
        NotDetected,                    // Peptide was not detected by MS2 or MBR. Only used in FlashLFQ Results when performing peptide-level quantification
    }
}
