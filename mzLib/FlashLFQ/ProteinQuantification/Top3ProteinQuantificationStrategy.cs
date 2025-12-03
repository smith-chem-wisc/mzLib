using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// Top3 protein quantification strategy: protein intensity is calculated as the sum of 
    /// the top 3 peptide intensities in each file
    /// </summary>
    public class Top3ProteinQuantificationStrategy<TPeptide, TProteinGroup, TSpectraFile> 
        : IProteinQuantificationStrategy<TPeptide, TProteinGroup, TSpectraFile>
        where TPeptide : IQuantifiablePeptide
        where TProteinGroup : IQuantifiableProteinGroup
        where TSpectraFile : IQuantifiableSpectraFile
    {
        private readonly int _topN;

        public Top3ProteinQuantificationStrategy(int topN = 3)
        {
            _topN = topN;
        }

        public void QuantifyProteins(
            Dictionary<TProteinGroup, List<TPeptide>> proteinGroupToPeptides,
            List<TSpectraFile> spectraFiles,
            bool useSharedPeptides)
        {
            foreach (var kvp in proteinGroupToPeptides)
            {
                TProteinGroup proteinGroup = kvp.Key;
                List<TPeptide> peptidesForThisProtein = kvp.Value;

                foreach (TSpectraFile file in spectraFiles)
                {
                    // Get top N peptides by intensity in this file
                    double proteinIntensity = peptidesForThisProtein
                        .Where(p => p.ProteinGroups.Count() == 1 || useSharedPeptides)
                        .Select(p => p.GetIntensity(file))
                        .OrderByDescending(intensity => intensity)
                        .Take(_topN)
                        .Sum();

                    proteinGroup.SetIntensity(file, proteinIntensity);
                }
            }
        }
    }
}
