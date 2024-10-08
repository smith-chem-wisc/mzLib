namespace Readers
{
    public enum Software
    {
        Unspecified,
        MassSpecFile,

        // Deconvolution Results
        FLASHDeconv, // files tested were outputted from OpenMs3.0.0
        TopFD,       // files tested were outputted from v1.6.2

        // Search Results
        MetaMorpheus,
        MaxQuant,
        Toppic,
        MsFragger, // files tested were from fragpipe v21.1
        MsPathFinderT,
        Crux,
        ProteomeDiscoverer,
        ProsightPD,
        Chimerys,

        // Quantification Results
        FlashLFQ,
    }
}
