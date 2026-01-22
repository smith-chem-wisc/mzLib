using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace Omics.BioPolymer
{
    /// <summary>
    /// Plain-language wrapper for a single VCF record (a line in a VCF file) with
    /// lightweight parsing of:
    /// - Reference and alternate allele strings
    /// - INFO (only passed through to <see cref="SnpEffAnnotation"/> for ANN-style annotations)
    /// - FORMAT column and per-sample genotype fields
    /// - Genotype (GT) tokens and Allelic Depth (AD) values
    /// - Simple zygosity classification per sample
    ///
    /// Design goals:
    /// - Fast, minimal allocation parsing for downstream proteomics / variant application.
    /// - Tolerant of missing data ('.') without throwing.
    /// - Avoids full VCF spec complexity (e.g., phased blocks, PL, GQ, allele remapping in multi-allelic normalization).
    ///
    /// Important assumptions / limitations:
    /// 1. The input line MUST be tab-delimited. Literal "\t" sequences will NOT be interpreted as tabs.
    /// 2. A valid VCF record is expected to contain at least the first 10 columns. If fewer are found, the constructor
    ///    returns early and most properties remain null / empty.
    /// 3. Only the ANN sub-field of INFO is parsed (via <see cref="SnpEffAnnotation"/>); all other INFO keys are ignored.
    /// 4. FORMAT fields are assumed to be consistent across all samples; mismatched token counts throw.
    /// 5. GT parsing:
    ///    - Splits on '/' or '|' and removes the separators.
    ///    - Missing alleles '.' are preserved in the parsed array.
    ///    - Unsupported allele indexes (>3) are still accepted if they appear (so long as they are numeric) – current validation allows 0–3 and '.'.
    /// 6. Zygosity rules:
    ///    - Only non-missing (not ".") allele symbols are considered.
    ///    - No called alleles ⇒ <see cref="Zygosity.Unknown"/>.
    ///    - One distinct called allele ⇒ <see cref="Zygosity.Homozygous"/>.
    ///    - More than one distinct called allele ⇒ <see cref="Zygosity.Heterozygous"/>.
    /// 7. Backward compatibility booleans (<see cref="Homozygous"/> / <see cref="Heterozygous"/>) are derived from the zygosity classification
    ///    and should be considered legacy conveniences. Prefer <see cref="ZygosityBySample"/>.
    ///
    /// Common usage pattern:
    /// <code>
    /// var vcf = new VariantCallFormat(vcLine);
    /// foreach (var (sampleId, gt) in vcf.Genotypes)
    /// {
    ///     var z = vcf.ZygosityBySample[sampleId];
    ///     var ad = vcf.AlleleDepths[sampleId];
    /// }
    /// </code>
    /// </summary>
    public class VariantCallFormat
    {
        /// <summary>
        /// Zygosity classification per sample, derived ONLY from called (non-missing) allele symbols.
        /// Missing-only genotype (e.g., "./.") ⇒ Unknown.
        /// </summary>
        public enum Zygosity { Unknown, Homozygous, Heterozygous }

        /// <summary>
        /// True when the provided line was truncated (< 10 VCF columns). In this case:
        /// - ReferenceAlleleString / AlternateAlleleString are null
        /// - AlleleIndex = -1
        /// - Info is a safe empty annotation (never null)
        /// - Format is an empty string
        /// - Genotypes / AlleleDepths / zygosity maps are empty
        /// </summary>
        public bool IsTruncated { get; }
        /// <summary>
        /// Original raw VCF line.
        /// </summary>
        public string Description { get; }

        /// <summary>
        /// REF allele text (may be null if constructor aborted).
        /// </summary>
        public string? ReferenceAlleleString { get; }

        /// <summary>
        /// ALT allele(s) comma-delimited (may be null if constructor aborted).
        /// </summary>
        public string? AlternateAlleleString { get; }

        /// <summary>
        /// Parsed snpEff-style annotation (ANN=*). All other INFO keys are ignored.
        /// </summary>
        public SnpEffAnnotation Info { get; }

        /// <summary>
        /// FORMAT column descriptor (e.g., "GT:AD:DP"). Used to parse sample columns.
        /// </summary>
        public string Format { get; }

        /// <summary>
        /// Per-sample genotype token arrays (GT split on '/' or '|').
        /// Keys are zero-based sample indices as strings ("0", "1", ...).
        /// </summary>
        public Dictionary<string, string[]> Genotypes { get; } = new();

        /// <summary>
        /// Per-sample AD (allele depth) string arrays (the raw comma-separated numeric tokens, excluding empty entries).
        /// Missing or invalid AD yields an empty array.
        /// </summary>
        public Dictionary<string, string[]> AlleleDepths { get; } = new();

        /// <summary>
        /// 1-based index of the allele referenced by ANN’s Allele (1..N for ALT, 0 for REF).
        /// -1 if the annotation's allele is missing or not found in ALT list.
        /// </summary>
        public int AlleleIndex { get; }

        /// <summary>
        /// Legacy: per-sample boolean flags indicating homozygosity.
        /// Prefer using <see cref="ZygosityBySample"/>.
        /// </summary>
        public Dictionary<string, bool> Homozygous { get; } = new();

        /// <summary>
        /// Legacy: per-sample boolean flags indicating heterozygosity.
        /// Prefer using <see cref="ZygosityBySample"/>.
        /// </summary>
        public Dictionary<string, bool> Heterozygous { get; } = new();

        /// <summary>
        /// Per-sample zygosity classification derived from non-missing genotype alleles.
        /// </summary>
        public Dictionary<string, Zygosity> ZygosityBySample { get; } = new();

        /// <summary>
        /// Construct from a single, tab-delimited VCF record.
        /// If fewer than 10 columns are present, parsing is aborted (object remains mostly unpopulated).
        /// </summary>
        /// <param name="description">Full raw VCF line (must contain actual tab characters).</param>
        public VariantCallFormat(string description)
        {
            if (description is null)
            {
                Description = string.Empty;
                ReferenceAlleleString = null;
                AlternateAlleleString = null;
                Info = new SnpEffAnnotation(string.Empty); // safe empty annotation
                Format = string.Empty;
                AlleleIndex = -1;
                IsTruncated = true;
                return;
            }

            Description = description;

            // Back-compat: if no real tabs are present but literal "\t" sequences are,
            // normalize them to actual tabs for parsing only. Leave Description intact.
            string parseLine = NormalizeTabsForParsing(description);

            // Parse description into VCF fields
            string[] vcfFields = parseLine.Split('\t');
            if (vcfFields.Length < 10)
            {
                ReferenceAlleleString = null;
                AlternateAlleleString = null;
                Info = new SnpEffAnnotation(string.Empty); // safe empty annotation
                Format = string.Empty;
                AlleleIndex = -1;
                IsTruncated = true;
                return;
            }
            ReferenceAlleleString = vcfFields[3];
            AlternateAlleleString = vcfFields[4];
            Info = new SnpEffAnnotation(vcfFields[7]);

            AlleleIndex = Info.Allele == null 
                ? -1
                : AlternateAlleleString.Split(',').ToList().IndexOf(Info.Allele) + 1; // returns 1-based index for ALT alleles, 0 if not found, -1 if Info.Allele is null

            // Format column tokens describe how to split each sample column
            Format = vcfFields[8];

            // Collect raw sample genotype strings (columns 9+)
            string[] genotypes = Enumerable
                .Range(9, vcfFields.Length - 9)
                .Select(i => vcfFields[i])
                .ToArray();

            // loop through genotypes for this variant (e.g. tumor and normal)
            // Parse each sample
            for (int individual = 0; individual < genotypes.Length; individual++)
            {
                var genotypeFields = GenotypeDictionary(Format.Trim(), genotypes[individual].Trim());

                // GT: split on '/' or '|' – separators removed intentionally.
                string[] gt = genotypeFields.TryGetValue("GT", out var gtString)
                    ? gtString.Split(new[] { '/', '|' })
                    : Array.Empty<string>();

                // Skip invalid or empty GT
                if (gt.Length == 0 || !GTvaluesAreValid(gt))
                {
                    continue;
                }

                // AD: optional – may be missing or contain '.' tokens
                int[] adDepths;
                string[] ad = genotypeFields.TryGetValue("AD", out var adString) && TryParseAD(adString, out adDepths)
                    ? adString.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
                    : Array.Empty<string>();

                string sampleKey = individual.ToString();
                Genotypes.Add(sampleKey, gt);
                AlleleDepths.Add(sampleKey, ad);

                // Zygosity classification: ignore '.' when counting distinct alleles
                var calledAlleles = gt.Where(a => a != ".").ToArray();
                Zygosity z;
                if (calledAlleles.Length == 0)
                {
                    z = Zygosity.Unknown;
                }
                else
                {
                    int distinctCalled = calledAlleles.Distinct().Count();
                    z = distinctCalled == 1 ? Zygosity.Homozygous : Zygosity.Heterozygous;
                }
                ZygosityBySample.Add(sampleKey, z);

                // Legacy boolean maps (retain for existing code paths)
                Homozygous.Add(sampleKey, z == Zygosity.Homozygous);
                Heterozygous.Add(sampleKey, z == Zygosity.Heterozygous);
            }
        }

        private static string NormalizeTabsForParsing(string line)
        {
            // Fast path: already contains real tabs
            if (line.IndexOf('\t') >= 0) return line;

            // Replace literal "\t" sequences for parsing only
            return line.IndexOf("\\t", StringComparison.Ordinal) >= 0
                ? line.Replace("\\t", "\t")
                : line;
        }

        /// <summary>
        /// Returns original string for the description
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return Description;
        }
        /// <summary>
        /// Equality is based solely on the original description string.
        /// </summary>
        public override bool Equals(object obj)
        {
            VariantCallFormat s = obj as VariantCallFormat;
            return s != null && s.Description == Description;
        }
        /// <summary>
        /// Hash code is derived from the original description (null-safe).
        /// </summary>
        public override int GetHashCode()
        {
            return (Description ?? "").GetHashCode();
        }

        /// <summary>
        /// Build a dictionary mapping FORMAT keys (e.g., GT, AD, DP) to the corresponding colon-delimited
        /// values from a single sample column. Throws if token counts differ.
        /// </summary>
        /// <param name="format">FORMAT column (e.g., "GT:AD:DP").</param>
        /// <param name="genotype">Sample column (e.g., "0/1:12,8:20").</param>
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

        /// <summary>
        /// Validate that all genotype tokens are drawn from the accepted set {0,1,2,3,.}.
        /// This is intentionally minimal; higher ALT indexes or symbolic alleles are not fully enforced here.
        /// </summary>
        public bool GTvaluesAreValid(string[] gt)
        {
            string[] validValues = { "0", "1", "2", "3", "." };
            return ValidationHelpers.TryValidateValues(gt.ToList(), validValues, out _);
        }

        /// <summary>
        /// Validate AD tokens: each must be "." or a non-negative integer.
        /// Empty AD arrays are considered invalid (if AD is present it should have content or '.').
        /// </summary>
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

        /// <summary>
        /// Attempt to parse AD into integer depths (excluding "." entries).
        /// Returns false if validation fails. On success, 'depths' contains only numeric values.
        /// </summary>
        public bool TryParseAD(string adString, out int[] depths)
        {
            depths = Array.Empty<int>();
            if (string.IsNullOrWhiteSpace(adString)) return false;

            var parts = adString.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries);
            if (!ADvaluesAreValid(parts)) return false;

            depths = parts.Where(p => p != ".").Select(int.Parse).ToArray();
            return true;
        }

        /// <summary>
        /// Shared validation helper for small, fixed vocabularies of acceptable string tokens.
        /// </summary>
        public static class ValidationHelpers
        {
            /// <summary>
            /// Returns true if all non-null, normalized values belong to the allowed set.
            /// Produces a distinct list of invalid tokens (if any).
            /// </summary>
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
                    seq
                        .Where(v => v is not null)
                        .Select(v => trim ? v!.Trim() : v!)
                        .Where(v => v.Length > 0);

                var normalized = Normalize(values);
                invalid = normalized
                    .Where(v => !allowed.Contains(v))
                    .Distinct(comparer)
                    .ToArray();
                return invalid.Length == 0;
            }
        }
    }
}