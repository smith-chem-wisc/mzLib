using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class ProteinWithAppliedVariants
        : Protein
    {
        public Protein Protein { get; }
        public List<SequenceVariation> AppliedSequenceVariations { get; }
        public string Individual { get; }

        public ProteinWithAppliedVariants(string variantBaseSequence, Protein protein, IEnumerable<SequenceVariation> appliedSequenceVariations, string individual)
            : base(variantBaseSequence, protein.Accession, organism: protein.Organism, gene_names: new List<Tuple<string, string>>(protein.GeneNames),
                  oneBasedModifications: protein.OneBasedPossibleLocalizedModifications.ToDictionary(x => x.Key, x => x.Value), 
                  proteolysisProducts: new List<ProteolysisProduct>(protein.ProteolysisProducts), name: protein.Name, full_name: protein.FullName, 
                  isDecoy: protein.IsDecoy, isContaminant: protein.IsContaminant, 
                  databaseReferences: new List<DatabaseReference>(protein.DatabaseReferences), sequenceVariations: new List<SequenceVariation>(protein.SequenceVariations),
                  disulfideBonds: new List<DisulfideBond>(protein.DisulfideBonds), databaseFilePath: protein.DatabaseFilePath)
        {
            Protein = protein;
            AppliedSequenceVariations = appliedSequenceVariations != null ? appliedSequenceVariations.ToList() : new List<SequenceVariation>();
            Individual = individual;
        }

        internal List<ProteinWithAppliedVariants> ApplyVariants(ProteinWithAppliedVariants protein, List<SequenceVariation> uniqueEffectsToApply)
        {
            // If there aren't any variants to apply, just return the base protein
            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<ProteinWithAppliedVariants> { protein };
            }

            List<ProteinWithAppliedVariants> proteins = new List<ProteinWithAppliedVariants>();
            foreach (SequenceVariation variant in uniqueEffectsToApply)
            {
                // Parse description into
                string[] vcfFields = variant.Description.Split('\t');
                if (vcfFields.Length < 10) { continue; }
                string referenceAlleleString = vcfFields[3];
                string alternateAlleleString = vcfFields[4];
                string info = vcfFields[7];
                string format = vcfFields[8];
                string[] genotypes = Enumerable.Range(9, vcfFields.Length - 9).Select(i => vcfFields[i]).ToArray();

                // loop through genotypes for this variant (e.g. tumor and normal)
                for (int i = 0; i < genotypes.Length; i++)
                {
                    if (Individual != null && Individual != i.ToString()) { continue; }
                    var genotypeFields = GenotypeDictionary(format.Trim(), genotypes[i].Trim());

                    // parse genotype
                    string[] gt = null;
                    if (genotypeFields.TryGetValue("GT", out string gtString)) { gt = gtString.Split('/'); }
                    if (gt == null) { continue; }

                    // parse allele depth (might be null, technically, but shouldn't be in most use cases)
                    string[] ad = null;
                    if (genotypeFields.TryGetValue("AD", out string adString)) { ad = adString.Split(','); }

                    // reference allele
                    if (gt.Contains("0"))
                    {
                        proteins.Add(new ProteinWithAppliedVariants(BaseSequence, Protein, AppliedSequenceVariations, i.ToString()));
                    }

                    // alternate allele
                    // TODO: recursively apply variants to create haplotypes and be wary of combinitorial explosion
                    if (!gt.All(x => x == "0"))
                    {
                        // check to see if there is incomplete indel overlap, which would lead to weird variant sequences
                        // complete overlap is okay, since it will be overwritten; this can happen if there are two alternate alleles,
                        //    e.g. reference sequence is wrong at that point
                        bool intersectsAppliedRegionIncompletely = AppliedSequenceVariations.Any(x => variant.Intersects(x) && !variant.Includes(x));
                        string seqBefore = BaseSequence.Substring(0, variant.OneBasedBeginPosition - 1);
                        string seqVariant = variant.VariantSequence;
                        int afterIdx = variant.OneBasedBeginPosition + variant.OriginalSequence.Length - 1;
                        if (intersectsAppliedRegionIncompletely)
                        {
                            // use original protein sequence for the remaining sequence
                            string seqAfter = Protein.BaseSequence.Length - afterIdx <= 0 ? "" : Protein.BaseSequence.Substring(afterIdx);
                            proteins.Add(new ProteinWithAppliedVariants(seqBefore + seqVariant + seqAfter, Protein, new[] { variant }, i.ToString()));
                        }
                        else
                        {
                            List<SequenceVariation> variations = AppliedSequenceVariations
                                .Where(x => !variant.Includes(x))
                                .Concat(new[] { variant })
                                .ToList();
                            // use the variant sequence for the remaining sequence
                            string seqAfter = BaseSequence.Length - afterIdx <= 0 ? "" : BaseSequence.Substring(afterIdx);
                            proteins.Add(new ProteinWithAppliedVariants(seqBefore + seqVariant + seqAfter, Protein, variations, i.ToString()));
                        }
                    }
                }
            }
            return proteins;
        }

        /// <summary>
        /// Gets a dictionary of the format (key) and fields (value) for a genotype
        /// </summary>
        /// <param name="format"></param>
        /// <param name="genotype"></param>
        /// <returns></returns>
        private static Dictionary<string, string> GenotypeDictionary(string format, string genotype)
        {
            Dictionary<string, string> genotypeDict = new Dictionary<string, string>();
            string[] formatSplit = format.Split(':');
            string[] genotypeSplit = genotype.Split(':');
            if (formatSplit.Length != genotypeSplit.Length)
            {
                throw new ArgumentException("Genotype format: " + format + " and genotype: " + genotype + " do not match -- they're not the same length");
            }
            return Enumerable.Range(0, formatSplit.Length).ToDictionary(x => formatSplit[x], x => genotypeSplit[x]);
        }
    }
}