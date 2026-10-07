using Omics.SequenceConversion;

namespace Readers.ProForma
{
    /// <summary>
    /// Registers ProForma 2.0 as a source and target format with a
    /// <see cref="SequenceConversionService"/>.
    /// <para>
    /// ProForma lives in <c>Readers</c> rather than <c>Omics</c> because it needs the
    /// TopDownProteomics SDK, which <c>Omics</c> does not reference; registration is therefore done at
    /// runtime, which <see cref="SequenceConversionService"/> supports by design. It is deliberately
    /// explicit rather than a module initializer: touching
    /// <see cref="SequenceConversionService.Default"/> loads the modification databases, and that cost
    /// should be paid by callers who convert sequences, not by every consumer of <c>Readers</c>.
    /// </para>
    /// <example>
    /// <code>
    /// ProFormaSequenceConversion.RegisterWithDefault();
    /// var proForma = SequenceConversionService.Default.Convert(
    ///     "PEP[Oxidation on M]TIDE", sourceFormat: "mzLib", targetFormat: "ProForma");
    /// </code>
    /// </example>
    /// </summary>
    public static class ProFormaSequenceConversion
    {
        private static readonly object RegistrationLock = new();
        private static bool _registeredWithDefault;

        /// <summary>
        /// Registers the ProForma parser, serializer, and the converters that pair ProForma with the
        /// formats already present on <paramref name="service"/>.
        /// </summary>
        public static void RegisterWith(SequenceConversionService service)
        {
            ArgumentNullException.ThrowIfNull(service);

            service.RegisterParser(ProFormaSequenceParser.Instance);
            service.RegisterSerializer(ProFormaSequenceSerializer.Instance);

            // ProForma -> ProForma is the canonicalizing round trip; the cross-format pairs let the
            // service convert without synthesizing a converter on demand.
            service.RegisterConverter(new SequenceConverter(ProFormaSequenceParser.Instance, ProFormaSequenceSerializer.Instance));
            service.RegisterConverter(new SequenceConverter(MzLibSequenceParser.Instance, ProFormaSequenceSerializer.Instance));
            service.RegisterConverter(new SequenceConverter(ProFormaSequenceParser.Instance, MzLibSequenceSerializer.Instance));
        }

        /// <summary>
        /// Registers ProForma with <see cref="SequenceConversionService.Default"/>. Idempotent, so it
        /// is safe to call from every entry point that needs ProForma conversion.
        /// </summary>
        public static void RegisterWithDefault()
        {
            if (_registeredWithDefault)
                return;

            lock (RegistrationLock)
            {
                if (_registeredWithDefault)
                    return;

                RegisterWith(SequenceConversionService.Default);
                _registeredWithDefault = true;
            }
        }
    }
}
