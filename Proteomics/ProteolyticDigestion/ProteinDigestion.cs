using System.Collections.Generic;
using System.Linq;
using Proteomics.Fragmentation;

namespace Proteomics.ProteolyticDigestion
{
    public class ProteinDigestion
    {
        /// <summary>
        /// Initializes digestion object
        /// </summary>
        /// <param name="digestionParams"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="variableModifications"></param>
        public ProteinDigestion(DigestionParams digestionParams, IEnumerable<Modification> allKnownFixedModifications, List<Modification> variableModifications)
        {
            DigestionParams = digestionParams;
            Protease = digestionParams.Protease;
            MaximumMissedCleavages = digestionParams.MaxMissedCleavages;
            InitiatorMethionineBehavior = digestionParams.InitiatorMethionineBehavior;
            MinPeptidesLength = digestionParams.MinPeptideLength;
            MaxPeptidesLength = digestionParams.MaxPeptideLength;
            AllKnownFixedModifications = allKnownFixedModifications ?? new List<Modification>();
            VariableModifications = variableModifications ?? new List<Modification>();
        }

        public Protease Protease { get; set; }
        public int MaximumMissedCleavages { get; set; }
        public DigestionParams DigestionParams { get; set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        public int MinPeptidesLength { get; set; }
        public int MaxPeptidesLength { get; set; }
        public IEnumerable<Modification> AllKnownFixedModifications { get; set; }
        public List<Modification> VariableModifications { get; set; }

        /// <summary>
        /// Gets peptides for semispecific digestion of a protein
        ///
        /// semi-specific search enters here...
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> SemiSpecificDigestion(Protein protein)
        {
            List<ProteolyticPeptide> intervals = new List<ProteolyticPeptide>();
            List<int> oneBasedIndicesToCleaveAfter = Protease.GetDigestionSiteIndices(protein.BaseSequence);

            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - MaximumMissedCleavages - 1; i++)
            {
                if (Protease.Retain(i, InitiatorMethionineBehavior, protein[0])
                    && Protease.OkayLength(oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], MinPeptidesLength, MaxPeptidesLength))
                {
                    intervals.Add(new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1],
                        oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], "semi"));
                }

                if (Protease.Cleave(i, InitiatorMethionineBehavior, protein[0])
                    && Protease.OkayLength(oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - 1, MinPeptidesLength, MaxPeptidesLength))
                {
                    intervals.Add(new ProteolyticPeptide(protein, 2, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1],
                        oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - 1, "semi:M cleaved"));
                }
            }

            int lastIndex = oneBasedIndicesToCleaveAfter.Count - 1;
            int maxIndex = MaximumMissedCleavages < lastIndex ? MaximumMissedCleavages : lastIndex;
            for (int i = 1; i <= maxIndex; i++)
            {
                if (DigestionParams.TerminusTypeSemiProtease == FragmentationTerminus.N) //tricky, it's N because we want the extra peptide at the C terminus |_
                {
                    if (Protease.OkayLength(oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], MinPeptidesLength, MaxPeptidesLength))
                    {
                        intervals.Add(new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[lastIndex - i] + 1, oneBasedIndicesToCleaveAfter[lastIndex],
                            oneBasedIndicesToCleaveAfter[lastIndex] - oneBasedIndicesToCleaveAfter[lastIndex - i], "semiN"));
                    }
                }
                else //TerminusType.C
                {
                    if (Protease.OkayLength(oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], MinPeptidesLength, MaxPeptidesLength))
                    {
                        intervals.Add(new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[0] + 1, oneBasedIndicesToCleaveAfter[i],
                            oneBasedIndicesToCleaveAfter[i] - oneBasedIndicesToCleaveAfter[0], "semiC"));
                    }
                }
            }

            // Also digest using the proteolysis product start/end indices
            intervals.AddRange(
                protein.ProteolysisProducts
                    .Where(proteolysisProduct => proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue &&
                    (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length))
                    .Select(proteolysisProduct => new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value,
                        0, proteolysisProduct.Type + " start")));

            return intervals.SelectMany(peptide => peptide.GetModifiedPeptides(AllKnownFixedModifications, DigestionParams, VariableModifications));
        }

        /// <summary>
        /// Gets peptides for specific protease digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> Digestion(Protein protein)
        {
            var intervals = Protease.GetDigestionIntervals(protein, MaximumMissedCleavages, InitiatorMethionineBehavior, MinPeptidesLength, MaxPeptidesLength);
            return intervals.SelectMany(peptide => peptide.GetModifiedPeptides(AllKnownFixedModifications, DigestionParams, VariableModifications));
        }
    }
}