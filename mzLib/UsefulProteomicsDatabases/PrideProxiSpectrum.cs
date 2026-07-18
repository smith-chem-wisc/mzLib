using System;
using System.Collections.Generic;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One spectrum returned by a PROXI query, keyed by a USI (Universal Spectrum Identifier). PROXI is the
    /// PSI standard REST API for retrieving a spectrum by USI; PRIDE serves it at
    /// <c>/pride/proxi/archive/v0.1/spectra?usi={usi}&amp;resultType=full</c>. A plain data object populated by
    /// JSON deserialization: the parallel <see cref="Mzs"/>/<see cref="Intensities"/> peak arrays convert to an
    /// mzLib <see cref="MassSpectrometry.MzSpectrum"/> via <see cref="PrideArchiveExtensions.ToMzSpectrum"/>, and
    /// the controlled-vocabulary <see cref="Attributes"/> (charge, precursor m/z, ms level, scan number,
    /// instrument, ...) are modeled with <see cref="CvParam"/>.
    /// </summary>
    /// <remarks>
    /// PROXI is repository-agnostic — the same shape is served by PeptideAtlas, MassIVE, and ProteomeCentral —
    /// so this DTO is not PRIDE-specific beyond where <see cref="PrideArchiveClient"/> points. The wire cvParams
    /// carry <c>accession</c>/<c>name</c>/(optional) <c>value</c> but no <c>cvLabel</c>; the accession prefix
    /// ("MS:", "PRIDE:", "UO:", ...) identifies the source ontology, so match on <see cref="CvParam.Accession"/>.
    /// </remarks>
    public class PrideProxiSpectrum
    {
        /// <summary>The spectrum's read status as reported by PROXI (e.g. "READABLE").</summary>
        public string Status { get; set; } = string.Empty;

        /// <summary>The Universal Spectrum Identifier this spectrum was resolved from.</summary>
        public string Usi { get; set; } = string.Empty;

        /// <summary>The m/z values of the spectrum's peaks. Parallel to <see cref="Intensities"/>.</summary>
        public double[] Mzs { get; set; } = Array.Empty<double>();

        /// <summary>The intensity values of the spectrum's peaks. Parallel to <see cref="Mzs"/>.</summary>
        public double[] Intensities { get; set; } = Array.Empty<double>();

        /// <summary>
        /// Spectrum and dataset metadata as controlled-vocabulary terms — e.g. charge state (MS:1000041),
        /// selected ion m/z (MS:1000744), ms level (MS:1000511), scan number (MS:1003057), instrument model.
        /// </summary>
        public List<CvParam> Attributes { get; set; } = new();
    }
}
