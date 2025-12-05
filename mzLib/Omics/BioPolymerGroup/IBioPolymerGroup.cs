using MassSpectrometry;
using Omics.Modifications;

namespace Omics.BioPolymerGroup
{
    internal interface IBioPolymerGroup : IEquatable<IBioPolymerGroup>
    {
        bool IsDecoy { get; }

        bool IsContaminant { get; }

        List<SpectraFileInfo> FilesForQuantification { get; set; }

        HashSet<IBioPolymer> BioPolymers { get; set; }

        string BioPolymerGroupName { get; }

        double BioPolymerGroupScore { get; set; }

        HashSet<IBioPolymerWithSetMods> AllBioPolymerWithSetMods { get; set; }

        HashSet<IBioPolymerWithSetMods> UniqueBioPolymerWithSetMods { get; set; }

        HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR { get; set; }

        double QValue { get; set; }

        double BestBiopolymerWithSetsModQValue { get; set; }

        double BestBioPolymerWithSetsModScore { get; set; }

        List<string> ModsInfo { get; }

        Dictionary<SpectraFileInfo, double> IntensitiesByFile { get; set; }

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
        void Score();
        void MergeWith(IBioPolymerGroup otherBioPolymerGroup);
        IBioPolymerGroup ConstructSubsetBioPolymerGroup(string fullFilePath, List<SilacLabel> silacLabels = null);
    }
}
