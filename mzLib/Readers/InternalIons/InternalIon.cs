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
        public double TicNormalizedIntensity { get; set; }
        public double TotalIonCurrent { get; set; }

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

        public char InternalNTerminalAA => string.IsNullOrEmpty(InternalSequence) ? '-' : InternalSequence[0];
        public char InternalCTerminalAA => string.IsNullOrEmpty(InternalSequence) ? '-' : InternalSequence[^1];

        public string FullModifiedSequence { get; set; } = string.Empty;
        public string ModificationsInInternalFragment { get; set; } = string.Empty;
        public bool HasModifiedResidue => CommonBiologicalModCount > 0 || CommonArtifactModCount > 0 ||
                                          MetalModCount > 0 || LessCommonModCount > 0;

        public bool PassesMassAccuracyFilter
        {
            get
            {
                if (double.IsNaN(MassErrorPpm)) return false;
                double abs = System.Math.Abs(MassErrorPpm);
                return (abs < 5.0) || (abs < 15.0 && NormalizedIntensity > 0.10);
            }
        }

        // ========== B/Y TERMINAL ION FEATURES ==========

        public double BIonIntensityAtNTerm { get; set; }
        public double YIonIntensityAtCTerm { get; set; }
        public bool HasMatchedBIonAtNTerm { get; set; }
        public bool HasMatchedYIonAtCTerm { get; set; }
        public double BYProductScore { get; set; }

        public int DistanceFromCTerm { get; set; }

        /// <summary>
        /// Length of the full peptide sequence.
        /// </summary>
        public int PeptideLength { get; set; }

        /// <summary>
        /// Relative position of the fragment's C-terminal end within the peptide.
        /// Computed as DistanceFromCTerm / PeptideLength.
        /// Range: 0.0 (fragment ends at peptide C-terminus) to ~1.0 (fragment ends near N-terminus).
        /// </summary>
        public double RelativeDistanceFromCTerm { get; set; }

        public double MaxTerminalIonIntensity { get; set; }
        public bool HasBothTerminalIons { get; set; }

        public int BasicResiduesInBIonSpan { get; set; }
        public int BasicResiduesInYIonSpan { get; set; }

        // ========== SCORER FEATURES ==========

        public bool IsProlineAtInternalNTerminus { get; set; }
        public bool IsTerminalRescue { get; set; }
        public double NTerminalFlankingHydrophobicity { get; set; }

        // ========== MODIFICATION CATEGORY FEATURES ==========

        public int CommonBiologicalModCount { get; set; }
        public int CommonArtifactModCount { get; set; }
        public int MetalModCount { get; set; }
        public int LessCommonModCount { get; set; }
        public bool HasPhosphorylation { get; set; }
        public bool HasMetalOnTerminalAcidic { get; set; }

        public static string[] GetHeaderNames() => new[]
        {
            nameof(PeptideSequence), nameof(InternalSequence), nameof(StartResidue), nameof(EndResidue),
            nameof(FragmentLength), nameof(TheoreticalMass), nameof(ObservedMass), nameof(MassError),
            nameof(MassErrorPpm), nameof(NormalizedIntensity), nameof(TicNormalizedIntensity), nameof(TotalIonCurrent),
            nameof(LocalIntensityRank), nameof(PrecursorCharge), nameof(CollisionEnergy), nameof(PeptidePEP),
            nameof(PeptideScore), nameof(HasProlineAtEitherTerminus), nameof(HasAspartateAtEitherTerminus),
            nameof(NTerminalFlankingResidue), nameof(CTerminalFlankingResidue), nameof(NumberOfBasicResidues),
            nameof(IsDecoy), nameof(SourceFile), nameof(ScanNumber), nameof(IsIsobaricAmbiguous),
            nameof(InternalNTerminalAA), nameof(InternalCTerminalAA), nameof(FullModifiedSequence),
            nameof(ModificationsInInternalFragment), nameof(HasModifiedResidue), nameof(PassesMassAccuracyFilter),
            nameof(BIonIntensityAtNTerm), nameof(YIonIntensityAtCTerm), nameof(HasMatchedBIonAtNTerm),
            nameof(HasMatchedYIonAtCTerm), nameof(BYProductScore), nameof(DistanceFromCTerm),
            nameof(PeptideLength), nameof(RelativeDistanceFromCTerm),
            nameof(MaxTerminalIonIntensity), nameof(HasBothTerminalIons),
            nameof(BasicResiduesInBIonSpan), nameof(BasicResiduesInYIonSpan),
            nameof(IsProlineAtInternalNTerminus), nameof(IsTerminalRescue), nameof(NTerminalFlankingHydrophobicity),
            nameof(CommonBiologicalModCount), nameof(CommonArtifactModCount), nameof(MetalModCount),
            nameof(LessCommonModCount), nameof(HasPhosphorylation), nameof(HasMetalOnTerminalAcidic)
        };

        public string[] GetValues() => new[]
        {
            PeptideSequence, InternalSequence, StartResidue.ToString(), EndResidue.ToString(),
            FragmentLength.ToString(), TheoreticalMass.ToString("G17", CultureInfo.InvariantCulture),
            ObservedMass.ToString("G17", CultureInfo.InvariantCulture),
            MassError.ToString("G17", CultureInfo.InvariantCulture),
            MassErrorPpm.ToString("G17", CultureInfo.InvariantCulture),
            NormalizedIntensity.ToString("G17", CultureInfo.InvariantCulture),
            TicNormalizedIntensity.ToString("G17", CultureInfo.InvariantCulture),
            TotalIonCurrent.ToString("G17", CultureInfo.InvariantCulture),
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
            PeptideLength.ToString(),
            RelativeDistanceFromCTerm.ToString("G17", CultureInfo.InvariantCulture),
            MaxTerminalIonIntensity.ToString("G17", CultureInfo.InvariantCulture),
            HasBothTerminalIons.ToString(),
            BasicResiduesInBIonSpan.ToString(), BasicResiduesInYIonSpan.ToString(),
            IsProlineAtInternalNTerminus.ToString(), IsTerminalRescue.ToString(),
            NTerminalFlankingHydrophobicity.ToString("G17", CultureInfo.InvariantCulture),
            CommonBiologicalModCount.ToString(), CommonArtifactModCount.ToString(),
            MetalModCount.ToString(), LessCommonModCount.ToString(),
            HasPhosphorylation.ToString(), HasMetalOnTerminalAcidic.ToString()
        };
    }
}