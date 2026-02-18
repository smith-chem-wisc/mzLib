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

        /// <summary>
        /// True if this ion shares a theoretical mass (to 4 decimal places) with another internal ion in the same PSM.
        /// </summary>
        public bool IsIsobaricAmbiguous { get; set; }

        /// <summary>
        /// First amino acid of the internal fragment sequence.
        /// </summary>
        public char InternalNTerminalAA =>
            string.IsNullOrEmpty(InternalSequence) ? '-' : InternalSequence[0];

        /// <summary>
        /// Last amino acid of the internal fragment sequence.
        /// </summary>
        public char InternalCTerminalAA =>
            string.IsNullOrEmpty(InternalSequence) ? '-' : InternalSequence[^1];

        /// <summary>
        /// The complete modified peptide sequence string from the PSM.
        /// </summary>
        public string FullModifiedSequence { get; set; } = string.Empty;

        /// <summary>
        /// Semicolon-delimited list of modifications whose positions fall within
        /// [StartResidue, EndResidue] inclusive.
        /// </summary>
        public string ModificationsInInternalFragment { get; set; } = string.Empty;

        /// <summary>
        /// True if ModificationsInInternalFragment is non-empty.
        /// </summary>
        public bool HasModifiedResidue => !string.IsNullOrEmpty(ModificationsInInternalFragment);

        /// <summary>
        /// Intensity-conditioned quality filter:
        /// True if (|MassErrorPpm| &lt; 5.0) OR (|MassErrorPpm| &lt; 15.0 AND NormalizedIntensity &gt; 0.10)
        /// </summary>
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
            nameof(PassesMassAccuracyFilter)
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
            PassesMassAccuracyFilter.ToString()
        };
    }
}