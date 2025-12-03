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
        bool IsDecoy { get; }

        bool IsContaminant { get; }

        //List<SpectraFileInfo> FilesForQuantification { get; set; }

        HashSet<IBioPolymer> BioPolymers { get; set; }

        string BioPolymerGroupName { get; }

        double BioPolymerGroupScore { get; set; }

        HashSet<IBioPolymerWithSetMods> AllBioPolymerWithSetMods { get; set; }

        HashSet<IBioPolymerWithSetMods> UniqueBioPolymerWithSetMods { get; set; }

        //HashSet<SpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        List<double> SequenceCoverageFraction { get; }

        List<string> SequenceCoverageDisplayList { get; }

        List<string> SequenceCoverageDisplayListWithMods { get; }

        List<string> FragmentSequenceCoverageDisplayList { get; }

        double QValue { get; set; }

        double BestBiopolymerWithSetModQValue { get; set; }

        double BestBioPolymerWithSetModScore { get; set; }

        int CumulativeTarget { get; set; }

        int CumulativeDecoy { get; set; }

        bool DisplayModsOnBioPolymerWithSetMods { get; set; }

        List<string> ModsInfo { get; }

        //Dictionary<SpectraFileInfo, double> IntensitiesByFile { get; set; }

        List<IBioPolymer> ListOfBioPolymersOrderedByAccession { get; }

        bool IEquatable<IBioPolymerGroup>.Equals(IBioPolymerGroup? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (other.GetType() != GetType()) return false;
            return BioPolymerGroupName == other.BioPolymerGroupName;
        }

        string GetTabSeparatedHeader();
        string ToString();
        double Score();
        //CalculateSequenceCoverage();
        //MergeWith(IBioPolymerGroup otherBioPolymerGroup);
        BioPolymerGroup ConstructSubsetProteinGroup(string fullFilePath, List<SilacLabel> silacLabels = null);
    }
}
