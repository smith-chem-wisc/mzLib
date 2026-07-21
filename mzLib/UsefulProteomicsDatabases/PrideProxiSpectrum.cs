using System;
using System.Collections.Generic;
using MassSpectrometry;
using MzLibUtil;
using Newtonsoft.Json;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One spectrum returned by a PROXI query, keyed by a USI (Universal Spectrum Identifier). PROXI is the
    /// PSI standard REST API for retrieving a spectrum by USI; PRIDE serves it at
    /// <c>/pride/proxi/archive/v0.1/spectra?usi={usi}&amp;resultType=full</c>.
    ///
    /// This <b>is</b> an mzLib <see cref="MzSpectrum"/> — the wire's parallel <c>mzs</c>/<c>intensities</c>
    /// arrays become the spectrum's <see cref="MzSpectrum.XArray"/>/<see cref="MzSpectrum.YArray"/>, so a
    /// fetched PROXI spectrum can be handed directly to anything in mzLib that takes an
    /// <see cref="MzSpectrum"/>, with no conversion step. On top of that it carries what an
    /// <see cref="MzSpectrum"/> cannot: the <see cref="Usi"/> it was resolved from, its PROXI
    /// <see cref="Status"/>, and the controlled-vocabulary <see cref="Attributes"/> (charge, precursor m/z,
    /// ms level, scan number, instrument, ...) modeled as <see cref="CvParam"/>. Those attributes are read
    /// into a full <see cref="MsDataScan"/> by <see cref="PrideArchiveExtensions.ToMsDataScan"/>.
    /// </summary>
    /// <remarks>
    /// PROXI is repository-agnostic — the same shape is served by PeptideAtlas, MassIVE, and ProteomeCentral —
    /// so this DTO is not PRIDE-specific beyond where <see cref="PrideArchiveClient"/> points. The wire cvParams
    /// carry <c>accession</c>/<c>name</c>/(optional) <c>value</c> but no <c>cvLabel</c>; the accession prefix
    /// ("MS:", "PRIDE:", "UO:", ...) identifies the source ontology, so match on <see cref="CvParam.Accession"/>.
    ///
    /// Because <see cref="MzSpectrum"/> has no parameterless constructor and keeps its peak arrays private to
    /// set, the peaks must arrive through the constructor rather than being assigned by the deserializer after
    /// the fact — hence the <see cref="JsonConstructorAttribute"/>. Newtonsoft matches the JSON keys to the
    /// constructor's parameter names case-insensitively, and supplies null for any key the server omits.
    /// </remarks>
    public class PrideProxiSpectrum : MzSpectrum
    {
        /// <summary>
        /// Builds a PROXI spectrum from its wire fields. Also the constructor Newtonsoft deserializes through.
        /// </summary>
        /// <param name="mzs">The m/z values of the spectrum's peaks. Null is treated as no peaks.</param>
        /// <param name="intensities">The intensities of the spectrum's peaks, parallel to <paramref name="mzs"/>.</param>
        /// <param name="usi">The Universal Spectrum Identifier this spectrum was resolved from.</param>
        /// <param name="status">The spectrum's read status as reported by PROXI (e.g. "READABLE").</param>
        /// <param name="attributes">The spectrum's controlled-vocabulary metadata terms.</param>
        /// <exception cref="MzLibException">The mzs and intensities arrays have different lengths.</exception>
        [JsonConstructor]
        public PrideProxiSpectrum(double[] mzs = null, double[] intensities = null, string usi = null,
            string status = null, List<CvParam> attributes = null)
            // shouldCopy: true because the sort below reorders the arrays in place; copying leaves a caller's
            // (or the deserializer's) arrays as they were, rather than reordering them behind its back.
            : base(ValidatedMzs(mzs, intensities, usi), intensities ?? Array.Empty<double>(), shouldCopy: true)
        {
            Usi = usi ?? string.Empty;
            Status = status ?? string.Empty;
            Attributes = attributes ?? new List<CvParam>();

            // MzSpectrum assumes ascending m/z and neither sorts nor validates, so sort defensively (PROXI
            // normally already returns them sorted). Sorting in place is safe here: the base constructor was
            // handed copies. Intensities ride along as the companion of the m/z key array.
            Array.Sort(XArray, YArray);
        }

        /// <summary>
        /// Coalesces the peak arrays and rejects a non-parallel pair, before the base constructor sees them.
        /// Static because it runs as part of the <c>base(...)</c> call, ahead of any instance state.
        /// </summary>
        private static double[] ValidatedMzs(double[] mzs, double[] intensities, string usi)
        {
            mzs ??= Array.Empty<double>();
            intensities ??= Array.Empty<double>();
            if (mzs.Length != intensities.Length)
                throw new MzLibException(
                    $"PROXI spectrum '{usi}' has {mzs.Length} m/z values but {intensities.Length} intensities; the peak arrays must be parallel.");
            return mzs;
        }

        /// <summary>The spectrum's read status as reported by PROXI (e.g. "READABLE").</summary>
        public string Status { get; }

        /// <summary>The Universal Spectrum Identifier this spectrum was resolved from.</summary>
        public string Usi { get; }

        /// <summary>
        /// Spectrum and dataset metadata as controlled-vocabulary terms — e.g. charge state (MS:1000041),
        /// selected ion m/z (MS:1000744), ms level (MS:1000511), scan number (MS:1003057), instrument model.
        /// A term may legitimately appear more than once (PRIDE sends one "instrument" term per instrument in
        /// the originating project), so match the first occurrence rather than expecting a single one.
        /// </summary>
        public List<CvParam> Attributes { get; }
    }
}
