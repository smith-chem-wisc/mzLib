using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Omics.SpectrumMatch
{
    public interface SpectrumMatchFromTsvHeader
    {
        // File and scan information
        const string FileName = "File Name";
        const string Ms2ScanNumber = "Scan Number";
        const string Ms2ScanRetentionTime = "Scan Retention Time";
        const string NumExperimentalPeaks = "Num Experimental Peaks";
        const string TotalIonCurrent = "Total Ion Current";
        const string PrecursorScanNum = "Precursor Scan Number";
        const string PrecursorCharge = "Precursor Charge";
        const string PrecursorMz = "Precursor MZ";
        const string PrecursorMass = "Precursor Mass";
        const string Score = "Score";
        const string DeltaScore = "Delta Score";
        const string Notch = "Notch";

        // Sequence information
        const string BaseSequence = "Base Sequence";
        const string FullSequence = "Full Sequence";
        const string EssentialSequence = "Essential Sequence";
        const string AmbiguityLevel = "Ambiguity Level";
        const string SpectrumMatchCount = "SM Count (unambiguous, <0.01 q-value)";
        const string Mods = "Mods";
        const string ModsChemicalFormulas = "Mods Chemical Formulas";
        const string ModsCombinedChemicalFormula = "Mods Combined Chemical Formula";
        const string NumVariableMods = "Num Variable Mods";
        const string MissedCleavages = "Missed Cleavages";
        const string MonoisotopicMass = "Monoisotopic Mass";
        const string MassDiffDa = "Mass Diff (Da)";
        const string MassDiffPpm = "Mass Diff (ppm)";
        const string Accession = "Accession";
        const string Name = "Name";
        const string GeneName = "Gene Name";
        const string OrganismName = "Organism Name";
        const string IntersectingSequenceVariations = "Intersecting Sequence Variations";
        const string IdentifiedSequenceVariations = "Identified Sequence Variations";
        const string SpliceSites = "Splice Sites";
        const string Contaminant = "Contaminant";
        const string Decoy = "Decoy";
        const string Description = "Description";
        const string StartAndEndResiduesInFullSequence = "Start and End Residues In Protein";
        const string PreviousMonomer = "Previous Monomer";
        const string NextMonomer = "Next Monomer";
        const string TheoreticalsSearched = "Theoreticals Searched";
        const string DecoyContaminantTarget = "Decoy/Contaminant/Target";
        const string MatchedIonSeries = "Matched Ion Series";
        const string MatchedIonMzRatios = "Matched Ion Mass-To-Charge Ratios";
        const string MatchedIonMassDiffDa = "Matched Ion Mass Diff (Da)";
        const string MatchedIonMassDiffPpm = "Matched Ion Mass Diff (Ppm)";
        const string MatchedIonIntensities = "Matched Ion Intensities";
        const string MatchedIonCounts = "Matched Ion Counts";

        // Scoring
        const string LocalizedScores = "Localized Scores";
        const string ImprovementPossible = "Improvement Possible";
        const string CumulativeTarget = "Cumulative Target";
        const string CumulativeDecoy = "Cumulative Decoy";
        const string CumulativeTargetNotch = "Cumulative Target Notch";
        const string CumulativeDecoyNotch = "Cumulative Decoy Notch";
        const string QValue = "QValue";
        const string QValueNotch = "QValue Notch";
        const string PEP = "PEP";
        const string PEP_QValue = "PEP_QValue";
        const string SpectralAngle = "Normalized Spectral Angle";
    }
}
