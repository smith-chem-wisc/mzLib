using Omics.SequenceConversion;

namespace Readers.ProForma
{
    /// <summary>
    /// Syntax description of HUPO-PSI ProForma 2.0 for the <see cref="SequenceConversionService"/>.
    /// <para>
    /// ProForma is not a bracket-tokenizable format in the sense the <c>SequenceParserBase</c> /
    /// <c>SequenceSerializerBase</c> helpers assume, so the ProForma parser and serializer implement
    /// <see cref="ISequenceParser"/> / <see cref="ISequenceSerializer"/> directly and delegate the
    /// string handling to Layer 1 (the TopDownProteomics reference parser/writer). This schema is
    /// therefore descriptive — it tells the framework what ProForma looks like; it is not the thing
    /// that parses it.
    /// </para>
    /// </summary>
    public sealed class ProFormaSequenceFormatSchema : SequenceFormatSchema
    {
        /// <summary>The format key this schema registers under.</summary>
        public const string ProFormaFormatName = "ProForma";

        public static ProFormaSequenceFormatSchema Instance { get; } = new();

        private ProFormaSequenceFormatSchema()
            : base(modOpen: '[', modClosed: ']', nTermSeparator: "-", cTermSeparator: "-")
        {
        }

        /// <inheritdoc />
        public override string FormatName => ProFormaFormatName;
    }
}
