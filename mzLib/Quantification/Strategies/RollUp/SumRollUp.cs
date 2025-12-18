using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;
using Quantification.Interfaces;

namespace Quantification.Strategies
{
    public class SumRollUp : RollUpStrategyBase, IRollUpStrategy
    {
        public override string Name => "Sum Roll-Up";
        // Implement roll-up methods here

        public SumRollUp() { }

        public override PeptideMatrix RollUpSpectralMatches(IExperimentalDesign experimentalDesign, List<ISpectralMatch> spectralMatches, List<IBioPolymerWithSetMods> peptides)
        {
            peptides.Sort((x, y) => x.FullSequence.CompareTo(y.FullSequence));
            var sampleInfos = GetOrderedSampleInfos(experimentalDesign, spectralMatches, out var filePathToArrayPositionDict);

            // Create empty matrix
            PeptideMatrix peptideMatrix = new PeptideMatrix(peptides, sampleInfos, experimentalDesign);

            // create empty array to store intensities as they are summed, before they're copied to the matrix
            double[] summedIntensities = new double[sampleInfos.Count];

            // Group spectral matches by peptide and then by file
            // The assumption when grouping by peptide is that no ambiguous spectral matches were passed in
            // If we start supporting ambiguous matches, we will need to adjust this logic
            foreach (var bioPolymerSpectralMatchGroup in spectralMatches
                .GroupBy(sm => sm.GetIdentifiedBioPolymersWithSetMods().First())
                .OrderBy(g => g.Key.FullSequence))
            {
                Array.Clear(summedIntensities, 0, summedIntensities.Length); // I'm pretty sure that values are copied when writing to the matrix,
                // so this is fine. If it doesn't work, we can create a new array each time instead.

                // Group by file
                var fileGroupedMatches = bioPolymerSpectralMatchGroup.GroupBy(g => g.FullFilePath);
                foreach (var fileGroup in fileGroupedMatches)
                {
                    // Each file is associated with one column (for LFQ) or multiple columns (for TMT) in the matrix
                    // The arrayPositions dict enables mapping between file paths and their corresponding array positions in the matrix
                    if (!filePathToArrayPositionDict.TryGetValue(fileGroup.Key, out var arrayPositions))
                    {
                        throw new KeyNotFoundException($"File path '{fileGroup.Key}' not found in file path to array position dictionary.");
                    }

                    // Finally, iterate through the spectral matches from the same file and sum their intensity
                    foreach (var sm in fileGroup)
                    {
                        // If no quant values are present, skip this SM
                        if (sm.QuantValues == null)
                        {
                            continue; 
                        }
                        
                        // Copy the quant values to the appropriate position in the summedIntensities array
                        for (var i = 0; i < arrayPositions.Count; i++)
                        {
                            var arrayPosition = arrayPositions[i];
                            summedIntensities[arrayPosition] += sm.QuantValues[i];
                        }
                    }
                }

                // After summing intensities for this peptide, write to the matrix
                peptideMatrix.SetRow(bioPolymerSpectralMatchGroup.Key, summedIntensities);
            }

            return peptideMatrix;
        }

        public override ProteinMatrix RollUpPeptides(PeptideMatrix peptides, List<IBioPolymerGroup> proteins)
        {
            throw new NotImplementedException();
            //ProteinMatrix result = new ProteinMatrix();
            //foreach (var peptide in peptides)
            //{


            //}
            //return result;
        }
    }
}
