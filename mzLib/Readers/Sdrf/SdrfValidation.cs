namespace Readers
{
    /// <summary>
    /// How badly a validation finding matters.
    ///
    /// The split is not cosmetic. A validator whose every rule is fatal is useless against real
    /// data: the 1,236-file curated corpus at bigbio/sdrf-annotated-datasets -- the closest thing
    /// the community has to a reference set -- violates parts of its own specification. A rule that
    /// fires on most curated files is a wrong rule, not 1,236 wrong files.
    /// </summary>
    internal enum SdrfValidationSeverity
    {
        /// <summary>
        /// The document deviates from the specification but can still be consumed and joined:
        /// a mixed-case column name, an empty cell where a reserved word belongs, a missing
        /// recommended column. Worth reporting, never worth refusing.
        /// </summary>
        Warning,

        /// <summary>
        /// The document cannot be reliably consumed: a row whose cells do not line up with the
        /// header, a missing column that everything keys on, two rows that are not distinguishable.
        /// </summary>
        Error
    }

    /// <summary>
    /// One validation finding, located as precisely as the rule allows.
    /// </summary>
    /// <param name="Severity">See <see cref="SdrfValidationSeverity"/>.</param>
    /// <param name="Rule">
    /// Stable identifier for the rule that fired (e.g. "MandatoryColumn", "RowKeyUniqueness").
    /// Stable so that callers can suppress a rule, and so that counting findings by rule across a
    /// corpus is meaningful.
    /// </param>
    /// <param name="Message">Human-readable description, including the offending value.</param>
    /// <param name="RowIndex">
    /// Zero-based index into <see cref="ResultFile{SdrfRow}.Results"/>, or null for a document-level
    /// finding. Note this is NOT the line number in the file: line number is RowIndex + 2, because
    /// the header occupies line 1 and file lines are 1-based.
    /// </param>
    /// <param name="ColumnName">The column involved, or null when the finding is not column-specific.</param>
    internal sealed record SdrfValidationMessage(
        SdrfValidationSeverity Severity,
        string Rule,
        string Message,
        int? RowIndex = null,
        string? ColumnName = null)
    {
        /// <summary>1-based line number in the source file, or null for a document-level finding.</summary>
        public int? LineNumber => RowIndex is null ? null : RowIndex + 2;

        public override string ToString()
        {
            string where = LineNumber is null ? "document" : $"line {LineNumber}";
            string column = ColumnName is null ? "" : $" [{ColumnName}]";
            return $"{Severity} {Rule} ({where}{column}): {Message}";
        }
    }

    /// <summary>
    /// The outcome of validating one SDRF document.
    /// </summary>
    internal sealed class SdrfValidationResult
    {
        public SdrfValidationResult(IEnumerable<SdrfValidationMessage> messages)
        {
            Messages = messages?.ToList() ?? throw new ArgumentNullException(nameof(messages));
        }

        public IReadOnlyList<SdrfValidationMessage> Messages { get; }

        public IEnumerable<SdrfValidationMessage> Errors =>
            Messages.Where(m => m.Severity == SdrfValidationSeverity.Error);

        public IEnumerable<SdrfValidationMessage> Warnings =>
            Messages.Where(m => m.Severity == SdrfValidationSeverity.Warning);

        /// <summary>True when nothing fatal was found. Warnings do not make a document invalid.</summary>
        public bool IsValid => !Errors.Any();

        public override string ToString() =>
            IsValid && Messages.Count == 0
                ? "valid"
                : $"{Errors.Count()} error(s), {Warnings.Count()} warning(s)";
    }
}
