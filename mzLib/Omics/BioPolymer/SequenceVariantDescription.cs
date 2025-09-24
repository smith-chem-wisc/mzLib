using System;
using System.Collections.Generic;
using System.Linq;

namespace Omics.BioPolymer
{
    public class SequenceVariantDescription
    {

        // Example VCF line with snpEff annotation:
        // 1   50000000   .   A   G   .   PASS   ANN=G||||||||||||||||   GT:AD:DP   1/1:30,30:30

        // --- VCF Standard Columns ---
        //
        // CHROM (1)      → Chromosome name (here, chromosome 1).
        // POS (50000000) → 1-based position of the variant (50,000,000).
        // ID (.)         → Variant identifier. "." means no ID (e.g., not in dbSNP).
        // REF (A)        → Reference allele in the reference genome (A).
        // ALT (G)        → Alternate allele observed in reads (G).
        // QUAL (.)       → Variant call quality score (Phred-scaled). "." means not provided.
        // FILTER (PASS)  → Indicates if the call passed filtering. "PASS" = high confidence.
        //
        // --- INFO Column ---
        //
        // INFO (ANN=...) holds snpEff annotation data.
        // ANN format is:
        //   Allele | Effect | Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID |
        //   Transcript_Biotype | Rank | HGVS.c | HGVS.p | cDNA_pos/cDNA_len |
        //   CDS_pos/CDS_len | AA_pos/AA_len | Distance | Errors/Warnings
        //
        // In this case: ANN=G||||||||||||||||
        //   - Allele = G
        //   - All other fields are empty → snpEff did not predict any functional impact
        //     (likely intergenic or unannotated region).
        //
        // --- FORMAT Column ---
        //
        // FORMAT (GT:AD:DP) defines how to read the sample column(s):
        //   GT → Genotype
        //   AD → Allele depth (number of reads supporting REF and ALT)
        //   DP → Read depth (total reads covering the site)
        //
        // --- SAMPLE Column ---
        //
        // Sample entry: 1/1:30,30:30
        //   GT = 1/1 → Homozygous ALT genotype (both alleles = G) ****************SEE NOTE BELOW******************
        //   AD = 30,30 → Read counts: REF=A has 30 reads, ALT=G has 30 reads
        //                (⚠ usually homozygous ALT would have few/no REF reads;
        //                 this may be caller-specific behavior or a quirk.)
        //   DP = 30   → Total coverage at this site = 30 reads
        //                (⚠ note AD sums to 60, which does not match DP.
        //                 This discrepancy is common in some callers.)
        //
        // --- Overall Summary ---
        // Variant at chr1:50000000 changes A → G.
        // The sample is homozygous for the ALT allele (G).
        // Variant passed filters, but no functional annotation from snpEff.

        // VCF GT (Genotype) Reference Key
        // --------------------------------
        //
        // Numbers correspond to alleles:
        //   0 = REF allele
        //   1 = first ALT allele
        //   2 = second ALT allele
        //   3 = third ALT allele (and so on)
        //
        // Symbols:
        //   / = unphased (we don't know which allele is on which chromosome)
        //   | = phased (we know which allele is on which haplotype)
        //   . = missing allele (no call)
        //
        // Common cases:
        // GT     Meaning                                   Example (REF=A, ALT=G)
        // 0/0    Homozygous reference                     A/A
        // 0/1    Heterozygous (REF + first ALT)           A/G
        // 1/0    Heterozygous (same as 0/1)               G/A
        // 1/1    Homozygous first ALT                     G/G
        // ././   Missing genotype                         -
        // 0|1    Phased heterozygous                      A on hap1, G on hap2
        // 1|0    Phased heterozygous (opposite phase)     G on hap1, A on hap2
        // .|1    One missing, one ALT                     missing/G
        // 0|.    One REF, one missing                     A/missing
        //
        // Multi-allelic examples (REF=A, ALT=G,T):
        // GT     Meaning                                   Example
        // 0/2    Heterozygous (REF + second ALT)          A/T
        // 1/2    Heterozygous (two different ALTs)        G/T
        // 2/2    Homozygous second ALT                    T/T
        // 0/3    Heterozygous (REF + third ALT)           A/[3rd ALT]
        // 2/3    Heterozygous (second + third ALT)        T/[3rd ALT]
        // 3/3    Homozygous third ALT                     [3rd ALT]/[3rd ALT]

        // VCF AD (Allelic Depths) and DP (Read Depth) Reference Key
        // ---------------------------------------------------------
        //
        // FORMAT field definitions:
        //   AD = Allelic depths for the ref and alt alleles in the order listed
        //   DP = Read depth (total number of reads covering the site)
        //
        // AD details:
        // - AD is usually represented as comma-separated integers.
        // - First value = reads supporting REF allele.
        // - Subsequent values = reads supporting each ALT allele in order.
        // - Example (REF=A, ALT=G):
        //     AD=35,12   -> 35 reads support A, 12 reads support G
        // - Example (REF=A, ALT=G,T):
        //     AD=40,5,10 -> 40 reads support A, 5 support G, 10 support T
        //
        // DP details:
        // - DP gives the total read depth across the site (may be equal to sum of AD, but not always).
        // - Sometimes DP includes low-quality or unfiltered reads that are not in AD.
        // - Example:
        //     AD=35,12, DP=47  -> total 47 reads, 35 REF, 12 ALT (0 reads mapped but not counted in AD)
        //     AD=40,5,10, DP=55 -> total 55 reads, 40 REF, 5 ALT1, 10 ALT2
        //
        // Special cases:
        // - AD=0,0 or DP=0  -> no reads cover this site.
        // - Missing values may be represented as "."
        //
        // Summary:
        //   AD helps you see how many reads support each allele individually.
        //   DP tells you the overall depth of coverage at the variant site.

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
            ReferenceAlleleString = vcfFields[3];
            AlternateAlleleString = vcfFields[4];
            Info = new SnpEffAnnotation(vcfFields[7]);
            AlleleIndex = Info.Allele == null ? -1 : AlternateAlleleString.Split(',').ToList().IndexOf(Info.Allele) + 1; // reference is zero
            Format = vcfFields[8];
            string[] genotypes = Enumerable.Range(9, vcfFields.Length - 9).Select(i => vcfFields[i]).ToArray();

            // loop through genotypes for this variant (e.g. tumor and normal)
            for (int individual = 0; individual < genotypes.Length; individual++)
            {
                var genotypeFields = GenotypeDictionary(Format.Trim(), genotypes[individual].Trim());

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
            }
        }

        public string Description { get; }
        public string ReferenceAlleleString { get; }
        public string AlternateAlleleString { get; }
        public SnpEffAnnotation Info { get; }
        public string Format { get; }
        public Dictionary<string, bool> Homozygous { get; } = new Dictionary<string, bool>();
        public Dictionary<string, bool> Heterozygous { get; } = new Dictionary<string, bool>();
        public Dictionary<string, string[]> Genotypes { get; } = new Dictionary<string, string[]>();
        public Dictionary<string, string[]> AlleleDepths { get; } = new Dictionary<string, string[]>();
        public int AlleleIndex { get; }

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
            return s != null && s.Description == Description;
        }

        public override int GetHashCode()
        {
            return (Description ?? "").GetHashCode();
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