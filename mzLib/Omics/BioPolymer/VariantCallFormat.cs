using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Omics.BioPolymer
{
    public class VariantCallFormat
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

        public VariantCallFormat(string description)
        {
            Description = description;

            // FIX: Split on actual tab characters instead of the literal sequence "\t"
            // Old (buggy): description.Split(new[] { @"\t" }, StringSplitOptions.None);
            string[] vcfFields = description.Split('\t');

            if (vcfFields.Length < 10)
            {
                ReferenceAlleleString = null;
                AlternateAlleleString = null;
                Info = null;
                Format = null;
                return;
            }

            ReferenceAlleleString = vcfFields[3];
            AlternateAlleleString = vcfFields[4];
            Info = new SnpEffAnnotation(vcfFields[7]);
            AlleleIndex = Info.Allele == null
                ? -1
                : AlternateAlleleString.Split(',').ToList().IndexOf(Info.Allele) + 1;
            Format = vcfFields[8];
            string[] genotypes = Enumerable.Range(9, vcfFields.Length - 9).Select(i => vcfFields[i]).ToArray();

            for (int individual = 0; individual < genotypes.Length; individual++)
            {
                var genotypeFields = GenotypeDictionary(Format.Trim(), genotypes[individual].Trim());

                string[] gt = genotypeFields.TryGetValue("GT", out var gtString)
                    ? gtString.Split(new[] { '/', '|' }, StringSplitOptions.RemoveEmptyEntries)
                    : Array.Empty<string>();

                if (gt.IsNullOrEmpty() && !GTvaluesAreValid(gt))
                {
                    continue;
                }

                int[] adDepths;
                string[] ad = genotypeFields.TryGetValue("AD", out var adString) && TryParseAD(adString, out adDepths)
                    ? adString.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                    : Array.Empty<string>();

                Genotypes.Add(individual.ToString(), gt);
                AlleleDepths.Add(individual.ToString(), ad);
                Homozygous.Add(individual.ToString(), gt.Distinct().Count() == 1);
                Heterozygous.Add(individual.ToString(), gt.Distinct().Count() > 1);
            }
        }

        public string Description { get; }
        public string? ReferenceAlleleString { get; }
        public string? AlternateAlleleString { get; }
        public SnpEffAnnotation Info { get; }
        public string Format { get; }
        public Dictionary<string, bool> Homozygous { get; } = new();
        public Dictionary<string, bool> Heterozygous { get; } = new();
        public Dictionary<string, string[]> Genotypes { get; } = new();
        public Dictionary<string, string[]> AlleleDepths { get; } = new();
        public int AlleleIndex { get; }

        public override string ToString() => Description;

        public override bool Equals(object obj)
        {
            var s = obj as VariantCallFormat;
            return s != null && s.Description == Description;
        }

        public override int GetHashCode() => (Description ?? "").GetHashCode();

        internal static Dictionary<string, string> GenotypeDictionary(string format, string genotype)
        {
            string[] formatSplit = format.Split(':');
            string[] genotypeSplit = genotype.Split(':');
            if (formatSplit.Length != genotypeSplit.Length)
            {
                throw new ArgumentException("Genotype format: " + format + " and genotype: " + genotype + " do not match -- they're not the same length");
            }
            return Enumerable.Range(0, formatSplit.Length).ToDictionary(x => formatSplit[x], x => genotypeSplit[x]);
        }

        public bool GTvaluesAreValid(string[] gt)
        {
            string[] validValues = { "0", "1", "2", "3", "." };
            return ValidationHelpers.TryValidateValues(gt.ToList(), validValues, out _);
        }

        public bool ADvaluesAreValid(string[] ad)
        {
            if (ad is null || ad.Length == 0) return false;
            foreach (var token in ad)
            {
                var s = token?.Trim();
                if (string.IsNullOrEmpty(s)) return false;
                if (s == ".") continue;
                if (!int.TryParse(s, out var n) || n < 0) return false;
            }
            return true;
        }

        public bool TryParseAD(string adString, out int[] depths)
        {
            depths = Array.Empty<int>();
            if (string.IsNullOrWhiteSpace(adString)) return false;

            var parts = adString.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
            if (!ADvaluesAreValid(parts)) return false;

            depths = parts.Where(p => p != ".").Select(int.Parse).ToArray();
            return true;
        }

        public static class ValidationHelpers
        {
            public static bool TryValidateValues(
                IEnumerable<string?> values,
                IEnumerable<string> allowedValues,
                out string[] invalid,
                bool ignoreCase = true,
                bool trim = true)
            {
                var comparer = ignoreCase ? StringComparer.OrdinalIgnoreCase : StringComparer.Ordinal;
                var allowed = new HashSet<string>(allowedValues, comparer);

                IEnumerable<string> Normalize(IEnumerable<string?> seq) =>
                    seq.Where(v => v is not null)
                       .Select(v => trim ? v!.Trim() : v!)
                       .Where(v => v.Length > 0);

                var normalized = Normalize(values);
                invalid = normalized.Where(v => !allowed.Contains(v))
                                    .Distinct(comparer)
                                    .ToArray();
                return invalid.Length == 0;
            }
        }
    }
}