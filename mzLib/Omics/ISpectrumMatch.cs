using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Omics.Fragmentation;

namespace Omics
{
    public interface ISpectrumMatch
    {
        string FullSequence { get; }
        int Ms2ScanNumber { get; }
        string FileNameWithoutExtension { get; }
        int PrecursorScanNum { get; }
        int PrecursorCharge { get; }
        double PrecursorMz { get; }
        double PrecursorMass { get; }
        double Score { get; }
        string ProteinAccession { get; }
        double? SpectralAngle { get; }
        List<MatchedFragmentIon> MatchedIons { get; }
        double QValue { get; }
        double PEP { get; }
        double PEP_QValue { get; }
        double? TotalIonCurrent { get; }
        double? DeltaScore { get; }
        string Notch { get; }
        string BaseSeq { get; }
        string EssentialSeq { get; }
        string AmbiguityLevel { get; }
        string MissedCleavage { get; }
        string PeptideMonoMass { get; }
        string MassDiffDa { get; }
        string MassDiffPpm { get; }
        string ProteinName { get; }
        string GeneName { get; }
        string OrganismName { get; }
        string IntersectingSequenceVariations { get; }
        string IdentifiedSequenceVariations { get; }
        //string PeptideDescription { get; }
        string StartAndEndResidues { get; }
        string PreviousMonomer { get; }
        string NextMonomer { get; }
        string DecoyContamTarget { get; }
        double? QValueNotch { get; }
        double? RetentionTime { get; }

        string ToString();
    }
}
