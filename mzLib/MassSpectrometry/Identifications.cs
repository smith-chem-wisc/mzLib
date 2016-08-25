using Spectra;

namespace MassSpectrometry
{
    public interface Identifications
    {
        int Count { get; }

        Tolerance parentTolerance { get; }

        Tolerance fragmentTolerance { get; }

        bool isDecoy(int matchIndex);

        int ms2spectrumIndex(int matchIndex);

        double calculatedMassToCharge(int matchIndex);

        double experimentalMassToCharge(int matchIndex);

        string PeptideSequenceWithoutModifications(int matchIndex);

        int chargeState(int matchIndex);

        int NumModifications(int matchIndex);

        int modificationLocation(int matchIndex, int i);

        string modificationDictionary(int matchIndex, int i);

        string modificationAcession(int matchIndex, int i);
    }
}
