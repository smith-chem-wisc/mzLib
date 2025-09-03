namespace Readers
{
    public enum Software
    {
        Unspecified,
        MassSpecFile,

        // Deconvolution
        FLASHDeconv, // files tested were outputted from OpenMs3.0.0
        TopFD,       // files tested were outputted from v1.6.2
        IsoDec,

        // Search
        MetaMorpheus,
        MaxQuant,
        Toppic,
        MsFragger, // files tested were from fragpipe v21.1
        MsPathFinderT,
        Crux,
        Dinosaur,
        ProteomeDiscoverer,
        ProsightPD,
        Chimerys,

        // Quantification
        FlashLFQ,
    }
}
