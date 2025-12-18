using Omics.Fragmentation;

namespace Omics.SpectralMatch;

/// <summary>
/// Interface for objects that can calculate and provide fragment-level sequence coverage.
/// Typically implemented by spectral match types that have matched fragment ion data
/// which can be used to determine which residues in the sequence are covered by fragments.
/// </summary>
public interface IHasSequenceCoverageFromFragments
{
    public string BaseSequence { get; }
    public List<MatchedFragmentIon> MatchedFragmentIons { get; protected set; }
    /// <summary>
    /// Positions in the biopolymer sequence (one-based) that are covered by fragment ions.
    /// Populated by <see cref="GetSequenceCoverage"/> when fragment coverage data is available.
    /// Returns null if coverage has not been calculated or no fragment data is available.
    /// </summary>
    public List<int> FragmentCoveragePositionInPeptide { get; protected set; }
    /// <summary>
    /// Calculates sequence coverage from fragment ions and populates 
    /// <see cref="FragmentCoveragePositionInPeptide"/> with one-based positions
    /// of residues that are covered by matched fragment ions.
    /// 
    /// Coverage is typically determined by identifying which residues have 
    /// fragments on both sides (N-terminal and C-terminal for peptides, 
    /// or 5' and 3' for nucleic acids).
    /// 
    /// Implementations without fragment ion data should leave the property null.
    /// This method may be called multiple times; implementations should 
    /// recalculate or return cached results as appropriate.
    /// </summary>
    public void GetSequenceCoverage();
}

public static class HasSequenceCoverageFromFragmentsExtensions
{
    /// <summary>
    /// Calculates the percentage of the sequence covered by fragment ions.
    /// Returns null if coverage has not been calculated or no fragment data is available.
    /// </summary>
    /// <param name="obj">The object implementing IHasSequenceCoverageFromFragments</param>
    /// <returns>Percentage of sequence covered by fragments, or null if not available</returns>
    public static List<int> GetSequenceCoverage(this IHasSequenceCoverageFromFragments obj)
    {
        if (string.IsNullOrEmpty(obj.BaseSequence) || !obj.MatchedFragmentIons.Any()) 
            return new List<int>();

        //Pull C terminal and N terminal Fragments and amino acid numbers
        var nTermFragmentAAPositions = obj.MatchedFragmentIons.Where(p =>
                p.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.N or FragmentationTerminus.FivePrime)
            .Select(j => j.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

        var cTermFragmentAAPositions = obj.MatchedFragmentIons.Where(p =>
                p.NeutralTheoreticalProduct.Terminus is FragmentationTerminus.C or FragmentationTerminus.ThreePrime)
            .Select(j => j.NeutralTheoreticalProduct.AminoAcidPosition).Distinct().ToList();

        //Create a hashset to store the covered amino acid positions
        HashSet<int> fragmentCoveredAminoAcids = new();

        //Check N term frags first
        if (nTermFragmentAAPositions.Any())
        {
            nTermFragmentAAPositions.Sort();

            //if the final NFragment is present, last AA is covered
            if (nTermFragmentAAPositions.Contains(obj.BaseSequence.Length - 1))
            {
                fragmentCoveredAminoAcids.Add(obj.BaseSequence.Length);
            }

            // if the first NFragment is present, first AA is covered
            if (nTermFragmentAAPositions.Contains(1))
            {
                fragmentCoveredAminoAcids.Add(1);
            }

            //Check all amino acids except for the last one in the list
            for (int i = 0; i < nTermFragmentAAPositions.Count - 1; i++)
            {
                //sequential AA, second one is covered
                if (nTermFragmentAAPositions[i + 1] - nTermFragmentAAPositions[i] == 1)
                {
                    fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1]);
                }

                //check to see if the position is covered from both directions, inclusive
                if (cTermFragmentAAPositions.Contains(nTermFragmentAAPositions[i + 1]))
                {
                    fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1]);
                }

                //check to see if the position is covered from both directions, exclusive
                if (cTermFragmentAAPositions.Contains(nTermFragmentAAPositions[i + 1] + 2))
                {
                    fragmentCoveredAminoAcids.Add(nTermFragmentAAPositions[i + 1] + 1);
                }
            }

        }

        //Check C term frags
        if (cTermFragmentAAPositions.Any())
        {
            cTermFragmentAAPositions.Sort();

            //if the second AA is present, the first AA is covered
            if (cTermFragmentAAPositions.Contains(2))
            {
                fragmentCoveredAminoAcids.Add(1);
            }

            //if the last AA is present, the final AA is covered
            if (cTermFragmentAAPositions.Contains(obj.BaseSequence.Length))
            {
                fragmentCoveredAminoAcids.Add(obj.BaseSequence.Length);
            }

            //check all amino acids except for the last one in the list
            for (int i = 0; i < cTermFragmentAAPositions.Count - 1; i++)
            {
                //sequential AA, the first one is covered
                if (cTermFragmentAAPositions[i + 1] - cTermFragmentAAPositions[i] == 1)
                {
                    fragmentCoveredAminoAcids.Add(cTermFragmentAAPositions[i]);
                }
            }
        }

        //store in PSM
        var fragmentCoveredAminoAcidsList = fragmentCoveredAminoAcids.ToList();
        fragmentCoveredAminoAcidsList.Sort();
        return fragmentCoveredAminoAcidsList;
    }
}

