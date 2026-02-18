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

        public bool IsIsobaricAmbiguous { get; set; }

        public char InternalNTerminalAA =>
            string.IsNullOrEmpty(InternalSequence) ? '-' : InternalSequence[0];

        public char InternalCTerminalAA =>
            string.IsNullOrEmpty(InternalSequence) ? '-' : InternalSequence[^1];

        public string FullModifiedSequence { get; set; } = string.Empty;
        public string ModificationsInInternalFragment { get; set; } = string.Empty;

        public bool HasModifiedResidue => !string.IsNullOrEmpty(ModificationsInInternalFragment);

        public bool PassesMassAccuracyFilter
        {
            get
            {
                if (double.IsNaN(MassErrorPpm))
                    return false;
                double absMassError = System.Math.Abs(MassErrorPpm);
                return (absMassError < 5.0) || (absMassError < 15.0 && NormalizedIntensity > 0.10);
            }
        }

        // ========== B/Y TERMINAL ION CORRELATION FEATURES ==========

        public double BIonIntensityAtNTerm { get; set; }
        public double YIonIntensityAtCTerm { get; set; }
        public bool HasMatchedBIonAtNTerm { get; set; }
        public bool HasMatchedYIonAtCTerm { get; set; }
        public double BYProductScore { get; set; }

        // ========== STEP 11 FEATURES ==========

        /// <summary>
        /// Distance from C-terminus: PeptideLength - EndResidue.
        /// 0 means the fragment ends at the last residue before the C-terminal cleavage site.
        /// </summary>
        public int DistanceFromCTerm { get; set; }

        /// <summary>
        /// Max(BIonIntensityAtNTerm, YIonIntensityAtCTerm).
        /// Captures the strongest terminal ion signal regardless of direction.
        /// </summary>
        public double MaxTerminalIonIntensity { get; set; }

        /// <summary>
        /// True if both b and y terminal ions are matched.
        /// </summary>
        public bool HasBothTerminalIons { get; set; }

        public static string[] GetHeaderNames() => new[]
        {
            nameof(PeptideSequence), nameof(InternalSequence), nameof(StartResidue), nameof(EndResidue),
            nameof(FragmentLength), nameof(TheoreticalMass), nameof(ObservedMass), nameof(MassError),
            nameof(MassErrorPpm), nameof(NormalizedIntensity), nameof(LocalIntensityRank), nameof(PrecursorCharge),
            nameof(CollisionEnergy), nameof(PeptidePEP), nameof(PeptideScore), nameof(HasProlineAtEitherTerminus),
            nameof(HasAspartateAtEitherTerminus), nameof(NTerminalFlankingResidue), nameof(CTerminalFlankingResidue),
            nameof(NumberOfBasicResidues), nameof(IsDecoy), nameof(SourceFile), nameof(ScanNumber),
            nameof(IsIsobaricAmbiguous), nameof(InternalNTerminalAA), nameof(InternalCTerminalAA),
            nameof(FullModifiedSequence), nameof(ModificationsInInternalFragment), nameof(HasModifiedResidue),
            nameof(PassesMassAccuracyFilter),
            nameof(BIonIntensityAtNTerm), nameof(YIonIntensityAtCTerm),
            nameof(HasMatchedBIonAtNTerm), nameof(HasMatchedYIonAtCTerm), nameof(BYProductScore),
            nameof(DistanceFromCTerm), nameof(MaxTerminalIonIntensity), nameof(HasBothTerminalIons)
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
            NumberOfBasicResidues.ToString(), IsDecoy.ToString(), SourceFile, ScanNumber,
            IsIsobaricAmbiguous.ToString(), InternalNTerminalAA.ToString(), InternalCTerminalAA.ToString(),
            FullModifiedSequence, ModificationsInInternalFragment, HasModifiedResidue.ToString(),
            PassesMassAccuracyFilter.ToString(),
            BIonIntensityAtNTerm.ToString("G17", CultureInfo.InvariantCulture),
            YIonIntensityAtCTerm.ToString("G17", CultureInfo.InvariantCulture),
            HasMatchedBIonAtNTerm.ToString(), HasMatchedYIonAtCTerm.ToString(),
            BYProductScore.ToString("G17", CultureInfo.InvariantCulture),
            DistanceFromCTerm.ToString(),
            MaxTerminalIonIntensity.ToString("G17", CultureInfo.InvariantCulture),
            HasBothTerminalIons.ToString()
        };
    }
}