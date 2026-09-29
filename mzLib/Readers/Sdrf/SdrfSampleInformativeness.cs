namespace Readers
{
    /// <summary>
    /// How much one SDRF says about its samples, as distinct from its assays.
    /// </summary>
    public enum SdrfSampleVerdict
    {
        /// <summary>
        /// Every check passes: a factor value varies, a sample characteristic other than organism is
        /// filled in, and the biological replicates are numbered apart. The file can drive a design.
        /// </summary>
        Informative,

        /// <summary>
        /// Some checks pass and some do not. Often legitimate -- a single-condition study has no
        /// factor to vary -- so what this means is the caller's decision, not this type's.
        /// </summary>
        Partial,

        /// <summary>
        /// No check passes. The sample half is reserved words, one repeated replicate number and no
        /// factor: a file generated from the data-file list, which for design purposes is the same
        /// as having no SDRF at all.
        /// </summary>
        Skeleton
    }

    /// <summary>
    /// The three sample checks on one document, with the column counts each was decided from.
    /// </summary>
    /// <param name="FactorValueColumns">Every <c>factor value[...]</c> column, not just the first.</param>
    /// <param name="SampleCharacteristicColumns">
    /// Every <c>characteristics[...]</c> column except organism and biological replicate. Organism is
    /// excluded because it is the one sample column a search can fill without a human, and so says
    /// nothing about whether anyone described the samples; biological replicate has its own check.
    /// </param>
    /// <param name="BiologicalReplicate">
    /// <c>characteristics[biological replicate]</c>, or null when the document lacks the column.
    /// </param>
    public sealed record SdrfSampleAssessment(
        IReadOnlyList<SdrfColumnCoverage> FactorValueColumns,
        IReadOnlyList<SdrfColumnCoverage> SampleCharacteristicColumns,
        SdrfColumnCoverage? BiologicalReplicate)
    {
        /// <summary>At least one factor value column holds two or more different real answers.</summary>
        public bool FactorValueVaries => FactorValueColumns.Any(c => c.DistinctValues >= 2);

        /// <summary>At least one sample characteristic, organism and replicate aside, holds a real answer.</summary>
        public bool SampleIsDescribed => SampleCharacteristicColumns.Any(c => c.Filled > 0);

        /// <summary>The biological replicate column holds two or more different real answers.</summary>
        public bool BiologicalReplicateVaries => BiologicalReplicate is { DistinctValues: >= 2 };

        public SdrfSampleVerdict Verdict =>
            (FactorValueVaries, SampleIsDescribed, BiologicalReplicateVaries) switch
            {
                (true, true, true) => SdrfSampleVerdict.Informative,
                (false, false, false) => SdrfSampleVerdict.Skeleton,
                _ => SdrfSampleVerdict.Partial
            };

        public override string ToString() =>
            $"{Verdict}: factor value varies {Mark(FactorValueVaries)}, sample described " +
            $"{Mark(SampleIsDescribed)}, biological replicate varies {Mark(BiologicalReplicateVaries)}";

        private static string Mark(bool passed) => passed ? "yes" : "no";
    }

    /// <summary>
    /// Decides, for ONE document, whether its sample half says anything -- the gate a pipeline needs
    /// before it trusts an SDRF to describe an experimental design.
    ///
    /// <see cref="SdrfCoverage"/> measures the same thing across a pooled corpus and cannot make this
    /// call on its own, for two reasons. It does not know sample columns from assay columns, and the
    /// assay columns are always filled, so a skeleton still looks half healthy. And its
    /// <see cref="SdrfColumnCoverage.IsUninformative"/> flags any column with one distinct value,
    /// which inside a single file is every organism column and every disease column of a
    /// single-condition study; gating on it would reject good files.
    ///
    /// So this asks three narrower questions, each over the counts <see cref="SdrfCoverage"/> already
    /// produces -- same reserved-word rule, same once-per-row counting of repeated columns:
    /// does a factor value vary, is any sample characteristic filled in, do the biological
    /// replicates differ. In bigbio/sdrf-annotated-datasets @ 4f823dcd, 194 of 1,236 curated files
    /// fail all three, and 193 of those carry no annotation-tool column, so the tool name cannot
    /// stand in for this check.
    ///
    /// Read-only, and it never throws on a bad document: the point is to describe one.
    /// </summary>
    public static class SdrfSampleInformativeness
    {
        private const string FactorValuePrefix = "factor value[";
        private const string CharacteristicsPrefix = "characteristics[";
        private const string Organism = "characteristics[organism]";
        private const string BiologicalReplicate = "characteristics[biological replicate]";

        public static SdrfSampleAssessment Assess(SdrfDocument document)
        {
            if (document is null) throw new ArgumentNullException(nameof(document));

            var coverage = SdrfCoverage.Measure(new SdrfCollection(new[] { document }, new[] { "document" }));

            // Prefixes are matched ignoring case: column names are otherwise compared ordinally, but
            // "Factor Value[" appears in two corpus files, and missing a factor column would turn an
            // informative file into a skeleton.
            var factors = coverage
                .Where(c => c.Column.StartsWith(FactorValuePrefix, StringComparison.OrdinalIgnoreCase))
                .OrderBy(c => c.Column, StringComparer.Ordinal)
                .ToList();

            var characteristics = coverage
                .Where(c => c.Column.StartsWith(CharacteristicsPrefix, StringComparison.OrdinalIgnoreCase)
                            && !string.Equals(c.Column, Organism, StringComparison.OrdinalIgnoreCase)
                            && !string.Equals(c.Column, BiologicalReplicate, StringComparison.OrdinalIgnoreCase))
                .OrderBy(c => c.Column, StringComparer.Ordinal)
                .ToList();

            var replicate = coverage.FirstOrDefault(c =>
                string.Equals(c.Column, BiologicalReplicate, StringComparison.OrdinalIgnoreCase));

            return new SdrfSampleAssessment(factors, characteristics, replicate);
        }
    }
}
