using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Omics.BioPolymerGroup
{
    internal interface IBioPolymerGroup : IEquatable<IBioPolymerGroup>
    {
        public bool IsDecoy { get; }

        public bool IsContaminant { get; }

        //public List<SpectraFileInfo> FilesForQuantification { get; set; }

        public HashSet<IBioPolymer> Biopolymers { get; set; }

        public string BiopolymerGroupName { get; }

        public double BiopolymerGroupScore { get; set; }

        public HashSet<IBioPolymerWithSetMods> AllPeptides { get; set; }

        public HashSet<IBioPolymerWithSetMods> UniquePeptides { get; set; }

        //public HashSet<SpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        public List<double> SequenceCoverageFraction { get; }

        public List<string> SequenceCoverageDisplayList { get; }

        public List<string> SequenceCoverageDisplayListWithMods { get; }

        public List<string> FragmentSequenceCoverageDisplayList { get; }

        public double QValue { get; set; }

        public double BestPeptideQValue { get; set; }

        public double BestPeptideScore { get; set; }

        public int CumulativeTarget { get; set; }

        public int CumulativeDecoy { get; set; }

        public bool DisplayModsOnPeptides { get; set; }

        public List<string> ModsInfo { get; }

        //public Dictionary<SpectraFileInfo, double> IntensitiesByFile { get; set; }

        public List<IBioPolymer> ListOfBiopolymersOrderedByAccession { get; }

        public string UniquePeptidesOutput { get; }
        public string SharedPeptidesOutput { get; }

        bool IEquatable<IBioPolymerGroup>.Equals(IBioPolymerGroup? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (other.GetType() != GetType()) return false;
            return BiopolymerGroupName == other.BiopolymerGroupName;
        }
    }
}
