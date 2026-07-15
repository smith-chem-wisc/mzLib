using System.Globalization;
using Omics.SequenceConversion;
using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// Serializes a <see cref="CanonicalSequence"/> to a ProForma 2.0 string, making ProForma a
    /// first-class target format for <see cref="SequenceConversionService"/>. Any format with a
    /// registered parser (mzLib, mass-shift, ProForma) can therefore be converted to standard
    /// ProForma notation without a bespoke converter.
    /// <para>
    /// Modifications are emitted as an ontology accession (<c>[UNIMOD:35]</c>) whenever one can be
    /// established — directly from the canonical UNIMOD id, or by resolving the modification through
    /// the configured <see cref="IModificationLookup"/> — and fall back to a name or mass descriptor
    /// otherwise. Term construction and string writing are delegated to
    /// <see cref="ProFormaConverter.BuildDescriptor"/> and Layer 1 (<see cref="ProFormaWriter"/>), so
    /// there is a single implementation of ProForma emission in mzLib.
    /// </para>
    /// </summary>
    public sealed class ProFormaSequenceSerializer : ISequenceSerializer
    {
        public static ProFormaSequenceSerializer Instance { get; } = new();

        /// <summary>
        /// Deferred so that merely registering this serializer does not force the modification
        /// databases to load.
        /// </summary>
        private readonly Lazy<IModificationLookup> _lookup;

        private ProFormaSequenceSerializer(IModificationLookup? lookup = null)
        {
            _lookup = lookup is null
                ? new Lazy<IModificationLookup>(() => MzLibModificationLookup.Instance)
                : new Lazy<IModificationLookup>(() => lookup);
        }

        /// <summary>Creates a serializer backed by a specific modification lookup.</summary>
        public static ProFormaSequenceSerializer WithLookup(IModificationLookup lookup)
        {
            ArgumentNullException.ThrowIfNull(lookup);
            return new ProFormaSequenceSerializer(lookup);
        }

        /// <inheritdoc />
        public string FormatName => ProFormaSequenceFormatSchema.ProFormaFormatName;

        /// <inheritdoc />
        public SequenceFormatSchema Schema => ProFormaSequenceFormatSchema.Instance;

        /// <inheritdoc />
        public IModificationLookup? ModificationLookup => _lookup.Value;

        /// <inheritdoc />
        public SequenceConversionHandlingMode HandlingMode => SequenceConversionHandlingMode.ThrowException;

        /// <summary>
        /// A <see cref="CanonicalSequence"/> is flat by construction, and ProForma can express every
        /// flat proteoform, so anything the canonical model can hold, ProForma can write.
        /// </summary>
        public bool CanSerialize(CanonicalSequence sequence) => !string.IsNullOrEmpty(sequence.BaseSequence);

        /// <summary>
        /// Resolve only when there is no accession to emit yet — an mzLib <c>Modification</c> carries
        /// the ontology references that <see cref="ProFormaConverter.BuildDescriptor"/> turns into a
        /// standard accession descriptor.
        /// </summary>
        public bool ShouldResolveMod(CanonicalModification mod) => !mod.IsResolved && !mod.UnimodId.HasValue;

        /// <inheritdoc />
        public string? Serialize(
            CanonicalSequence sequence,
            ConversionWarnings? warnings = null,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
        {
            warnings ??= new ConversionWarnings();

            if (!CanSerialize(sequence))
            {
                return SequenceConversionHelpers.HandleSerializerError(warnings, mode,
                    ConversionFailureReason.InvalidSequence, "Cannot write ProForma for an empty base sequence.");
            }

            List<Tdp.ProFormaDescriptor>? nTerm = null;
            List<Tdp.ProFormaDescriptor>? cTerm = null;
            var tags = new List<Tdp.ProFormaTag>();

            foreach (var mod in sequence.Modifications)
            {
                var enriched = ShouldResolveMod(mod) ? ModificationLookup?.TryResolve(mod) ?? mod : mod;
                var descriptor = BuildDescriptor(enriched);

                if (descriptor is null)
                {
                    var message = $"Modification '{mod.OriginalRepresentation}' at {mod.PositionType} carries no " +
                                  "accession, name, or mass, so it cannot be written as ProForma.";
                    warnings.AddIncompatibleItem(message);

                    if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
                    {
                        warnings.AddWarning($"Dropped {message}");
                        continue;
                    }

                    return SequenceConversionHelpers.HandleSerializerError(warnings, mode,
                        ConversionFailureReason.IncompatibleModifications, message);
                }

                switch (enriched.PositionType)
                {
                    case ModificationPositionType.NTerminus:
                        (nTerm ??= []).Add(descriptor);
                        break;
                    case ModificationPositionType.CTerminus:
                        (cTerm ??= []).Add(descriptor);
                        break;
                    default:
                        tags.Add(new Tdp.ProFormaTag(enriched.ResidueIndex!.Value, new[] { descriptor }));
                        break;
                }
            }

            tags.Sort((a, b) => a.ZeroBasedStartIndex.CompareTo(b.ZeroBasedStartIndex));

            var term = new Tdp.ProFormaTerm(sequence.BaseSequence, tags.Count > 0 ? tags : null, nTerm, cTerm,
                labileDescriptors: null, unlocalizedTags: null, tagGroups: null, globalModifications: null);

            return ProFormaWriter.Write(term);
        }

        /// <summary>
        /// Prefers a resolved mzLib modification (which yields a real ontology accession), then a
        /// UNIMOD id, then a mass, then the raw name. Returns null when the modification carries none
        /// of these.
        /// </summary>
        private static Tdp.ProFormaDescriptor? BuildDescriptor(CanonicalModification mod)
        {
            if (mod.MzLibModification is not null)
                return ProFormaConverter.BuildDescriptor(mod.MzLibModification);

            if (mod.UnimodId.HasValue)
                return new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Identifier, Tdp.ProFormaEvidenceType.Unimod,
                    $"UNIMOD:{mod.UnimodId.Value}");

            if (!string.IsNullOrWhiteSpace(mod.OriginalRepresentation))
                return new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Name, mod.OriginalRepresentation);

            var mass = mod.EffectiveMass;
            if (mass.HasValue)
                return new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Mass,
                    mass.Value.ToString("+0.####;-0.####", CultureInfo.InvariantCulture));

            return null;
        }
    }
}
