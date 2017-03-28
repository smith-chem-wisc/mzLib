using MzLibUtil;

namespace MassSpectrometry
{
    public interface IIdentifications
    {
        #region Public Properties

        int Count { get; }

        Tolerance ParentTolerance { get; }

        Tolerance FragmentTolerance { get; }

        #endregion Public Properties

        #region Public Methods

        string Ms2SpectrumID(int matchIndex);

        int ChargeState(int matchIndex);

        float[] MatchedIons(int matchIndex, int i);

        int MatchedIonCounts(int matchIndex, int i);

        string ProteinAccession(int matchIndex);

        string ProteinFullName(int matchIndex);

        int StartResidueInProtein(int matchIndex);

        int EndResidueInProtein(int matchIndex);

        bool IsDecoy(int matchIndex);

        bool PassThreshold(int matchIndex);

        double CalculatedMassToCharge(int matchIndex);

        double ExperimentalMassToCharge(int matchIndex);

        string PeptideSequenceWithoutModifications(int matchIndex);

        int NumModifications(int matchIndex);

        int ModificationLocation(int matchIndex, int i);

        string ModificationDictionary(int matchIndex, int i);

        string ModificationAcession(int matchIndex, int i);

        #endregion Public Methods
    }
}