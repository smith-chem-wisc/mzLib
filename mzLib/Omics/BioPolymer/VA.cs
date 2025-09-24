

namespace Omics.BioPolymer
{
    public static class VA
    {
        /// <summary>
        /// Creates a list of IBioPolymers of the same type as the original protein, each with applied variants from this protein.
        /// </summary>

        public static List<TBioPolymerType> GetVariantBioPolymers<TBioPolymerType>(this TBioPolymerType protein, int maxSequenceVariantsPerIsoform = 4, int minAlleleDepth = 1, int maxSequenceVariantIsoforms = 1)
            where TBioPolymerType : IHasSequenceVariants
        {
            List<TBioPolymerType> allBioplymers = new List<TBioPolymerType>() { protein};

            if (protein.SequenceVariations.All(v => v.AreValid()) && protein.SequenceVariations.Any(v => v.Description == null || v.Description.Genotypes.Count == 0))
            {
                // this is a protein with either no VCF lines or a mix of VCF and non-VCF lines
                allBioplymers.AddRange(ApplyAllVariantCombinations(protein, protein.SequenceVariations, maxSequenceVariantsPerIsoform, maxSequenceVariantIsoforms).ToList());
            }
            return allBioplymers;
        }
        /// <summary>
        /// Applies all possible combinations of the provided SequenceVariation list to the base TBioPolymerType object,
        /// starting with the fewest single variations and up to the specified maximum number of combinations.
        /// </summary>

        public static IEnumerable<TBioPolymerType> ApplyAllVariantCombinations<TBioPolymerType>(
            TBioPolymerType baseBioPolymer,
            List<SequenceVariation> variations,
            int maxSequenceVariantsPerIsoform,
            int maxSequenceVariantIsoforms)
            where TBioPolymerType : IHasSequenceVariants
        {
            int count = 0;

            // Always yield the base biopolymer first
            yield return baseBioPolymer;
            count++;
            //if (count >= maxSequenceVariantsPerIsoform)
            //    yield break;

            int n = variations.Count;
            // generate combinations of isoforms but limit the number of variants per isoform
            for (int size = 1; size <= maxSequenceVariantsPerIsoform; size++)
            {
                foreach (var combo in GetCombinations(variations, size))
                {
                    // break if we've reached the maximum number of isoforms
                    if (count >= maxSequenceVariantIsoforms)
                        yield break;
                    if (!ValidCombination(combo.ToList()))
                        continue;
                    var result = baseBioPolymer;
                    foreach (var variant in combo)
                    {
                        result = ApplySingleVariant(variant, result, string.Empty);
                    }
                    if (result != null)
                    {
                        yield return result;
                        count++;

                    }
                }
            }
        }
        /// <summary>
        /// Generates all possible combinations of the specified size from the input list.
        /// </summary>
        /// <param name="variations">List of SequenceVariation objects to combine. Assumed not null or empty.</param>
        /// <param name="size">The size of each combination.</param>
        /// <returns>
        /// An IEnumerable of IList&lt;SequenceVariation&gt; representing each combination.
        /// </returns>
        private static IEnumerable<IList<SequenceVariation>> GetCombinations(List<SequenceVariation> variations, int size)
        {
            int n = variations.Count;
            var indices = new int[size];
            for (int i = 0; i < size; i++) indices[i] = i;

            while (true)
            {
                var combo = new List<SequenceVariation>(size);
                for (int i = 0; i < size; i++)
                    combo.Add(variations[indices[i]]);
                yield return combo;

                int pos = size - 1;
                while (pos >= 0 && indices[pos] == n - size + pos)
                    pos--;
                if (pos < 0) break;
                indices[pos]++;
                for (int i = pos + 1; i < size; i++)
                    indices[i] = indices[i - 1] + 1;
            }
        }
    }
}
