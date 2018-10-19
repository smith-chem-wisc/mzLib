using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteomics
{
    public class SequenceVariantDescription
    {
        public SequenceVariantDescription(string description)
        {
            Description = description;
            if (description == null)
            {
                return;
            }

            // Parse description into
            string[] vcfFields = description.Split(new[] { @"\t" }, StringSplitOptions.None);
            if (vcfFields.Length < 10) { return; }
            string referenceAlleleString = vcfFields[3];
            string alternateAlleleString = vcfFields[4];
            string info = vcfFields[7];
            string format = vcfFields[8];
            string[] genotypes = Enumerable.Range(9, vcfFields.Length - 9).Select(i => vcfFields[i]).ToArray();

            // loop through genotypes for this variant (e.g. tumor and normal)
            for (int individual = 0; individual < genotypes.Length; individual++)
            {
                var genotypeFields = GenotypeDictionary(format.Trim(), genotypes[individual].Trim());

                // parse genotype
                string[] gt = null;
                if (genotypeFields.TryGetValue("GT", out string gtString)) { gt = gtString.Split('/'); }
                if (gt == null) { continue; }

                // parse allele depth (might be null, technically, but shouldn't be in most use cases)
                string[] ad = null;
                if (genotypeFields.TryGetValue("AD", out string adString)) { ad = adString.Split(','); }

                Genotypes.Add(individual.ToString(), gt);
                AlleleDepths.Add(individual.ToString(), ad);
                Homozygous.Add(individual.ToString(), gt.Distinct().Count() == 1);
                Heterozygous.Add(individual.ToString(), gt.Distinct().Count() > 1);
                GenotypeAlleleDepthMap.Add(individual.ToString(), Enumerable.Range(0, gt.Length).Select(x => (gt[x], ad[x])).ToArray());
            }
        }

        public string Description { get; }
        public string ReferenceAlleleString { get; }
        public string AlternateAlleleString { get; }
        public string Info { get; }
        public string Format { get; }
        public Dictionary<string, bool> Homozygous { get; } = new Dictionary<string, bool>();
        public Dictionary<string, bool> Heterozygous { get; } = new Dictionary<string, bool>();
        public Dictionary<string, string[]> Genotypes { get; } = new Dictionary<string, string[]>();
        public Dictionary<string, string[]> AlleleDepths { get; } = new Dictionary<string, string[]>();
        public Dictionary<string, (string, string)[]> GenotypeAlleleDepthMap { get; } = new Dictionary<string, (string, string)[]>();

        /// <summary>
        /// Returns original string for the description
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return Description;
        }

        public override bool Equals(object obj)
        {
            SequenceVariantDescription s = obj as SequenceVariantDescription;
            return s != null
                && (s.Description == null && Description == null || s.Description.Equals(Description));
        }

        public override int GetHashCode()
        {
            return Description.GetHashCode();
        }

        /// <summary>
        /// Gets a dictionary of the format (key) and fields (value) for a genotype
        /// </summary>
        /// <param name="format"></param>
        /// <param name="genotype"></param>
        /// <returns></returns>
        internal static Dictionary<string, string> GenotypeDictionary(string format, string genotype)
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