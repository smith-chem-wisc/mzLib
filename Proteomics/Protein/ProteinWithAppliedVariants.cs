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

        public ProteinWithAppliedVariants(string variantBaseSequence, Protein protein, IEnumerable<SequenceVariation> appliedSequenceVariations,
            IEnumerable<ProteolysisProduct> applicableProteolysisProducts, IDictionary<int, List<Modification>> oneBasedModifications, string individual)
            : base(variantBaseSequence,
                  protein.Accession + (appliedSequenceVariations == null ? "" : "_" + CombineSimpleStrings(appliedSequenceVariations)),
                  organism: protein.Organism,
                  geneNames: new List<Tuple<string, string>>(protein.GeneNames),
                  oneBasedModifications: oneBasedModifications.ToDictionary(x => x.Key, x => x.Value),
                  proteolysisProducts: new List<ProteolysisProduct>(applicableProteolysisProducts),
                  name: protein.Name + (appliedSequenceVariations == null ? "" : " variant:" + CombineDescriptions(appliedSequenceVariations)),
                  fullName: protein.FullName + (appliedSequenceVariations == null ? "" : " variant:" + CombineDescriptions(appliedSequenceVariations)),
                  isDecoy: protein.IsDecoy,
                  isContaminant: protein.IsContaminant,
                  databaseReferences: new List<DatabaseReference>(protein.DatabaseReferences),
                  sequenceVariations: new List<SequenceVariation>(protein.SequenceVariations),
                  disulfideBonds: new List<DisulfideBond>(protein.DisulfideBonds),
                  databaseFilePath: protein.DatabaseFilePath)
        {
            Protein = protein;
            AppliedSequenceVariations = appliedSequenceVariations != null ? appliedSequenceVariations.ToList() : new List<SequenceVariation>();
            Individual = individual;
        }

        /// <summary>
        /// Applies variant changes to protein sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="uniqueEffectsToApply"></param>
        /// <returns></returns>
        internal List<ProteinWithAppliedVariants> ApplyVariants(ProteinWithAppliedVariants protein, List<SequenceVariation> uniqueEffectsToApply)
        {
            // If there aren't any variants to apply, just return the base protein
            if (uniqueEffectsToApply.Count == 0)
            {
                return new List<ProteinWithAppliedVariants> { protein };
            }

            bool referenceAdded = false;
            List<ProteinWithAppliedVariants> proteins = new List<ProteinWithAppliedVariants>();
            foreach (SequenceVariation variant in uniqueEffectsToApply)
            {
                // Parse description into
                string[] vcfFields = variant.Description.Split(new[] { @"\t" }, StringSplitOptions.None);
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
                    if (gt.Contains("0") && !referenceAdded)
                    {
                        proteins.Add(new ProteinWithAppliedVariants(BaseSequence, Protein, AppliedSequenceVariations, ProteolysisProducts, OneBasedPossibleLocalizedModifications, i.ToString()));
                        referenceAdded = true; // only add the reference allele once
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
                        List<ProteolysisProduct> adjustedProteolysisProducts = AdjustProteolysisProductIndices(variant, ProteolysisProducts);
                        Dictionary<int, List<Modification>> adjustedModifications = AdjustModificationIndices(variant, OneBasedPossibleLocalizedModifications);
                        int afterIdx = variant.OneBasedBeginPosition + variant.OriginalSequence.Length - 1;
                        if (intersectsAppliedRegionIncompletely)
                        {
                            // use original protein sequence for the remaining sequence
                            string seqAfter = Protein.BaseSequence.Length - afterIdx <= 0 ? "" : Protein.BaseSequence.Substring(afterIdx);
                            proteins.Add(new ProteinWithAppliedVariants(seqBefore + seqVariant + seqAfter, Protein, new[] { variant }, adjustedProteolysisProducts, adjustedModifications, i.ToString()));
                        }
                        else
                        {
                            List<SequenceVariation> variations = AppliedSequenceVariations
                                .Where(x => !variant.Includes(x))
                                .Concat(new[] { variant })
                                .ToList();
                            // use this variant protein sequence for the remaining sequence
                            string seqAfter = BaseSequence.Length - afterIdx <= 0 ? "" : BaseSequence.Substring(afterIdx);
                            proteins.Add(new ProteinWithAppliedVariants(seqBefore + seqVariant + seqAfter, Protein, variations, adjustedProteolysisProducts, adjustedModifications, i.ToString()));
                        }
                    }
                }
            }
            return proteins;
        }

        /// <summary>
        /// Eliminates proteolysis products that overlap sequence variations. 
        /// Since frameshift indels are written across the remaining sequence, 
        /// this eliminates proteolysis products that conflict with large deletions and other structural variations.
        /// </summary>
        /// <param name="variants"></param>
        /// <param name="proteolysisProducts"></param>
        /// <returns></returns>
        internal List<ProteolysisProduct> AdjustProteolysisProductIndices(SequenceVariation variant, IEnumerable<ProteolysisProduct> proteolysisProducts)
        {
            List<ProteolysisProduct> products = new List<ProteolysisProduct>();
            if (proteolysisProducts == null) { return products; }
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;
            foreach (ProteolysisProduct p in proteolysisProducts.Where(p => p.OneBasedEndPosition.HasValue && p.OneBasedBeginPosition.HasValue))
            {
                // proteolysis product is entirely before the variant
                if (variant.OneBasedBeginPosition > p.OneBasedEndPosition)
                {
                    products.Add(p);
                }
                // proteolysis product straddles the variant, but the cleavage site(s) are still intact; the ends aren't considered cleavage sites
                else if ((p.OneBasedBeginPosition < variant.OneBasedBeginPosition || p.OneBasedBeginPosition == 1 || p.OneBasedBeginPosition == 2)
                    && (p.OneBasedEndPosition > variant.OneBasedEndPosition || p.OneBasedEndPosition == Protein.BaseSequence.Length)
                    && p.OneBasedEndPosition + sequenceLengthChange <= BaseSequence.Length)
                {
                    products.Add(new ProteolysisProduct(p.OneBasedBeginPosition, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                }
                // proteolysis product is after the variant
                else if (p.OneBasedBeginPosition > variant.OneBasedEndPosition
                    && p.OneBasedBeginPosition + sequenceLengthChange <= BaseSequence.Length && p.OneBasedEndPosition + sequenceLengthChange <= BaseSequence.Length)
                {
                    products.Add(new ProteolysisProduct(p.OneBasedBeginPosition + sequenceLengthChange, p.OneBasedEndPosition + sequenceLengthChange, p.Type));
                }
                else // sequence variant conflicts with proteolysis cleavage site (cleavage site was lost)
                {
                    continue;
                }
            }
            return products;
        }

        /// <summary>
        /// Adjusts modification indices.
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="modificationDictionary"></param>
        /// <returns></returns>
        internal Dictionary<int, List<Modification>> AdjustModificationIndices(SequenceVariation variant, IDictionary<int, List<Modification>> modificationDictionary)
        {
            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            if (modificationDictionary == null) { return mods; }
            int sequenceLengthChange = variant.VariantSequence.Length - variant.OriginalSequence.Length;
            foreach (KeyValuePair<int, List<Modification>> kv in modificationDictionary)
            {
                if (variant.OneBasedBeginPosition > kv.Key)
                {
                    mods.Add(kv.Key, kv.Value);
                }
                else if (variant.OneBasedEndPosition < kv.Key && kv.Key + sequenceLengthChange <= BaseSequence.Length)
                {
                    mods.Add(kv.Key + sequenceLengthChange, kv.Value);
                }
                else // sequence variant conflicts with modification site (modification site substitution)
                {
                    continue;
                }
            }
            return mods;
        }

        /// <summary>
        /// Format string to append to accession
        /// </summary>
        /// <param name="variations"></param>
        /// <returns></returns>
        internal static string CombineSimpleStrings(IEnumerable<SequenceVariation> variations)
        {
            return variations == null ? "" : string.Join("_", variations.Select(v => v.SimpleString()));
        }

        /// <summary>
        /// Format string to append to protein names
        /// </summary>
        /// <param name="variations"></param>
        /// <returns></returns>
        internal static string CombineDescriptions(IEnumerable<SequenceVariation> variations)
        {
            return variations == null ? "" : string.Join(", variant:", variations.Select(d => d.Description));
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