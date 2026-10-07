using System.Globalization;
using Omics.SequenceConversion;
using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// Parses a ProForma 2.0 string into a <see cref="CanonicalSequence"/>, making ProForma a
    /// first-class source format for <see cref="SequenceConversionService"/>.
    /// <para>
    /// The string handling is delegated to Layer 1 (<see cref="ProFormaReader"/>), which is lossless.
    /// The lossy step is Layer 1 -&gt; canonical: <see cref="CanonicalSequence"/> is a flat model
    /// (one modification per position, no groups/ranges/delocalization), so ProForma constructs that
    /// cannot be expressed in it — tag groups, position ranges, unlocalized, labile, and global
    /// modifications, and sequence ambiguity — are reported as incompatible rather than silently
    /// dropped. Use Layer 1 directly when losslessness matters.
    /// </para>
    /// </summary>
    public sealed class ProFormaSequenceParser : ISequenceParser
    {
        public static ProFormaSequenceParser Instance { get; } = new();

        private ProFormaSequenceParser() { }

        /// <inheritdoc />
        public string FormatName => ProFormaSequenceFormatSchema.ProFormaFormatName;

        /// <inheritdoc />
        public SequenceFormatSchema Schema => ProFormaSequenceFormatSchema.Instance;

        /// <summary>
        /// Fast heuristic used by <c>ParseAutoDetect</c>. Deliberately conservative: a bare sequence
        /// and mzLib's own <c>"[Type:Name on X]"</c> syntax are both left to the mzLib parser, so this
        /// only claims strings carrying a ProForma-distinctive token that also parse and fit the flat
        /// canonical model.
        /// </summary>
        public bool CanParse(string input)
        {
            if (string.IsNullOrWhiteSpace(input) || !LooksLikeProForma(input))
                return false;

            try
            {
                return DescribeUnrepresentable(ProFormaReader.Read(input)) is null;
            }
            catch (Tdp.ProFormaParseException)
            {
                return false;
            }
        }

        /// <inheritdoc />
        public CanonicalSequence? Parse(
            string input,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
        {
            warnings ??= new ConversionWarnings();

            if (string.IsNullOrWhiteSpace(input))
            {
                return SequenceConversionHelpers.HandleParserError(warnings, mode,
                    ConversionFailureReason.InvalidSequence, "ProForma input is null or empty.");
            }

            Tdp.ProFormaTerm term;
            try
            {
                term = ProFormaReader.Read(input);
            }
            catch (Tdp.ProFormaParseException ex)
            {
                return SequenceConversionHelpers.HandleParserError(warnings, mode,
                    ConversionFailureReason.InvalidSequence,
                    $"'{input}' is not valid ProForma 2.0: {ex.Message}");
            }

            var unrepresentable = DescribeUnrepresentable(term);
            if (unrepresentable is not null)
            {
                warnings.AddIncompatibleItem(unrepresentable);
                return SequenceConversionHelpers.HandleParserError(warnings, mode,
                    ConversionFailureReason.IncompatibleModifications,
                    $"ProForma term '{input}' cannot be represented as a CanonicalSequence: {unrepresentable}. " +
                    "Use the ProForma reader (Layer 1) directly to work with this term losslessly.");
            }

            var builder = new CanonicalSequenceBuilder(term.Sequence)
                .WithSourceFormat(FormatName);

            if (term.NTerminalDescriptors is { Count: > 0 })
            {
                var (representation, unimodId, mass) = Describe(term.NTerminalDescriptors, warnings, term);
                builder.AddNTerminalModification(representation, mass, formula: null, unimodId: unimodId);
            }

            if (term.Tags is not null)
            {
                foreach (var tag in term.Tags)
                {
                    var (representation, unimodId, mass) = Describe(tag.Descriptors, warnings, term);
                    builder.AddResidueModification(tag.ZeroBasedStartIndex, representation, mass,
                        formula: null, unimodId: unimodId);
                }
            }

            if (term.CTerminalDescriptors is { Count: > 0 })
            {
                var (representation, unimodId, mass) = Describe(term.CTerminalDescriptors, warnings, term);
                builder.AddCTerminalModification(representation, mass, formula: null, unimodId: unimodId);
            }

            // The flat model allows only one modification per position; ProForma does not.
            var errors = builder.Validate();
            if (errors.Count > 0)
            {
                foreach (var error in errors)
                    warnings.AddIncompatibleItem(error);

                return SequenceConversionHelpers.HandleParserError(warnings, mode,
                    ConversionFailureReason.IncompatibleModifications,
                    $"ProForma term '{input}' does not fit the canonical model: {string.Join(" ", errors)}");
            }

            return builder.Build();
        }

        /// <summary>
        /// Names the first ProForma construct in <paramref name="term"/> that the flat
        /// <see cref="CanonicalSequence"/> model cannot hold, or null when the term is representable.
        /// This is the Layer-1 → canonical boundary, stated explicitly so callers fail loud.
        /// </summary>
        internal static string? DescribeUnrepresentable(Tdp.ProFormaTerm term)
        {
            if (term.TagGroups is { Count: > 0 })
                return "tag groups (crosslinks / localization groups)";
            if (term.LabileDescriptors is { Count: > 0 })
                return "labile modifications";
            if (term.UnlocalizedTags is { Count: > 0 })
                return "unlocalized modifications";
            if (term.GlobalModifications is { Count: > 0 })
                return "global modifications";

            if (term.Tags is not null)
            {
                foreach (var tag in term.Tags)
                {
                    if (tag.ZeroBasedStartIndex != tag.ZeroBasedEndIndex)
                        return "position ranges";
                    if (tag.HasAmbiguousSequence)
                        return "sequence ambiguity";
                }
            }

            return null;
        }

        /// <summary>
        /// Collapses a ProForma descriptor list onto the single representation a
        /// <see cref="CanonicalModification"/> can hold, preferring an ontology accession (most
        /// resolvable) over a name over a mass. ProForma permits several descriptors on one tag
        /// (e.g. <c>[Phospho|UNIMOD:21]</c>); the canonical model does not, so the surplus is
        /// reported as a warning — the information is still intact in the Layer-1 term.
        /// </summary>
        private static (string Representation, int? UnimodId, double? Mass) Describe(
            IList<Tdp.ProFormaDescriptor> descriptors, ConversionWarnings warnings, Tdp.ProFormaTerm term)
        {
            Tdp.ProFormaDescriptor? preferred = null;
            int? unimodId = null;
            double? mass = null;

            foreach (var descriptor in descriptors)
            {
                if (descriptor.Key == Tdp.ProFormaKey.Identifier)
                {
                    preferred ??= descriptor;
                    unimodId ??= TryGetUnimodId(descriptor.Value);
                }
                else if (descriptor.Key == Tdp.ProFormaKey.Mass
                         && double.TryParse(descriptor.Value, NumberStyles.Float, CultureInfo.InvariantCulture, out var parsed))
                {
                    mass ??= parsed;
                }
            }

            preferred ??= descriptors.FirstOrDefault(d => d.Key == Tdp.ProFormaKey.Name)
                          ?? descriptors.FirstOrDefault();

            if (descriptors.Count > 1)
            {
                warnings.AddWarning(
                    $"ProForma tag [{string.Join("|", descriptors.Select(d => d.Value))}] in '{term.Sequence}' carries " +
                    $"{descriptors.Count} descriptors; the canonical model keeps one ('{preferred?.Value}').");
            }

            return (preferred?.Value ?? string.Empty, unimodId, mass);
        }

        /// <summary>Extracts 35 from "UNIMOD:35". Returns null for any other ontology.</summary>
        private static int? TryGetUnimodId(string? value)
        {
            const string prefix = "UNIMOD:";
            if (value is null || !value.StartsWith(prefix, StringComparison.OrdinalIgnoreCase))
                return null;

            return int.TryParse(value.AsSpan(prefix.Length), NumberStyles.Integer, CultureInfo.InvariantCulture, out var id)
                ? id
                : null;
        }

        /// <summary>
        /// True when the string carries a token only ProForma uses. Keeps auto-detection from
        /// claiming plain sequences or mzLib's "Name on X" motif syntax.
        /// </summary>
        private static bool LooksLikeProForma(string input)
        {
            if (input.Contains(" on ", StringComparison.OrdinalIgnoreCase))
                return false;

            if (input.IndexOfAny(['{', '<', '#']) >= 0)
                return true;

            foreach (var ontology in Ontologies)
            {
                if (input.Contains(ontology, StringComparison.OrdinalIgnoreCase))
                    return true;
            }

            // A bare mass tag, e.g. PEP[+15.9949]TIDE.
            for (int i = 0; i < input.Length - 2; i++)
            {
                if (input[i] == '[' && (input[i + 1] == '+' || input[i + 1] == '-') && char.IsDigit(input[i + 2]))
                    return true;
            }

            return false;
        }

        private static readonly string[] Ontologies =
            ["UNIMOD:", "MOD:", "RESID:", "XLMOD:", "GNO:", "Obs:", "Info:", "Formula:", "Glycan:"];
    }
}
