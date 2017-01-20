using Spectra;

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

        bool IsDecoy(int matchIndex);

        bool PassThreshold(int matchIndex);

        string Ms2SpectrumID(int matchIndex);

        double CalculatedMassToCharge(int matchIndex);

        double ExperimentalMassToCharge(int matchIndex);

        string PeptideSequenceWithoutModifications(int matchIndex);

        int ChargeState(int matchIndex);

        int NumModifications(int matchIndex);

        int ModificationLocation(int matchIndex, int i);

        string ModificationDictionary(int matchIndex, int i);

        string ModificationAcession(int matchIndex, int i);

        #endregion Public Methods

    }
}