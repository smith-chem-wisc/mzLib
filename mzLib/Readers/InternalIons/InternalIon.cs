using System.Globalization;
using System.Linq;

namespace Readers.InternalIons
{
    /// <summary>
    /// Represents a single annotated internal fragment ion extracted from a PSM.
    /// </summary>
    public class InternalFragmentIon
    {
        public string PeptideSequence { get; set; } = string.Empty;
        public string InternalSequence { get; set; } = string.Empty;
        public int StartResidue { get; set; }
        public int EndResidue { get; set; }
        public int FragmentLength => EndResidue - StartResidue + 1;
        public double TheoreticalMass { get; set; }
        public double ObservedMass { get; set; }
        public double MassError => ObservedMass - TheoreticalMass;
        public double MassErrorPpm => TheoreticalMass != 0 ? (MassError / TheoreticalMass) * 1e6 : double.NaN;
        public double NormalizedIntensity { get; set; }
        public double LocalIntensityRank { get; set; }
        public int PrecursorCharge { get; set; }
        public double CollisionEnergy { get; set; } = double.NaN;
        public double PeptidePEP { get; set; }
        public double PeptideScore { get; set; }

        public bool HasProlineAtEitherTerminus =>
            !string.IsNullOrEmpty(InternalSequence) &&
            (InternalSequence[0] == 'P' || InternalSequence[^1] == 'P');

        public bool HasAspartateAtEitherTerminus =>
            !string.IsNullOrEmpty(InternalSequence) &&
            (InternalSequence[0] == 'D' || InternalSequence[^1] == 'D');

        public char NTerminalFlankingResidue { get; set; } = '-';
        public char CTerminalFlankingResidue { get; set; } = '-';

        public int NumberOfBasicResidues =>
            string.IsNullOrEmpty(InternalSequence) ? 0 :
            InternalSequence.Count(c => c == 'K' || c == 'R' || c == 'H');

        public bool IsDecoy { get; set; }
        public string SourceFile { get; set; } = string.Empty;
        public string ScanNumber { get; set; } = string.Empty;

        public static string[] GetHeaderNames() => new[]
        {
            nameof(PeptideSequence), nameof(InternalSequence), nameof(StartResidue), nameof(EndResidue),
            nameof(FragmentLength), nameof(TheoreticalMass), nameof(ObservedMass), nameof(MassError),
            nameof(MassErrorPpm), nameof(NormalizedIntensity), nameof(LocalIntensityRank), nameof(PrecursorCharge),
            nameof(CollisionEnergy), nameof(PeptidePEP), nameof(PeptideScore), nameof(HasProlineAtEitherTerminus),
            nameof(HasAspartateAtEitherTerminus), nameof(NTerminalFlankingResidue), nameof(CTerminalFlankingResidue),
            nameof(NumberOfBasicResidues), nameof(IsDecoy), nameof(SourceFile), nameof(ScanNumber)
        };

        public string[] GetValues() => new[]
        {
            PeptideSequence, InternalSequence, StartResidue.ToString(), EndResidue.ToString(),
            FragmentLength.ToString(), TheoreticalMass.ToString("G17", CultureInfo.InvariantCulture),
            ObservedMass.ToString("G17", CultureInfo.InvariantCulture),
            MassError.ToString("G17", CultureInfo.InvariantCulture),
            MassErrorPpm.ToString("G17", CultureInfo.InvariantCulture),
            NormalizedIntensity.ToString("G17", CultureInfo.InvariantCulture),
            LocalIntensityRank.ToString("G17", CultureInfo.InvariantCulture),
            PrecursorCharge.ToString(),
            double.IsNaN(CollisionEnergy) ? "NaN" : CollisionEnergy.ToString("G17", CultureInfo.InvariantCulture),
            PeptidePEP.ToString("G17", CultureInfo.InvariantCulture),
            PeptideScore.ToString("G17", CultureInfo.InvariantCulture),
            HasProlineAtEitherTerminus.ToString(), HasAspartateAtEitherTerminus.ToString(),
            NTerminalFlankingResidue.ToString(), CTerminalFlankingResidue.ToString(),
            NumberOfBasicResidues.ToString(), IsDecoy.ToString(), SourceFile, ScanNumber
        };
    }
}