using System.Collections.Generic;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One sample characteristic of a PRIDE Archive project — a controlled-vocabulary key paired
    /// with one or more controlled-vocabulary values — as returned by the PRIDE Archive REST API
    /// (v3 <c>projects/{accession}</c>). A plain data object populated by JSON deserialization.
    /// </summary>
    /// <remarks>
    /// A single key carries a list of values because one characteristic can have several: the
    /// "organism" key (OBI:0100026) on a metaproteomics project lists every species in the sample.
    /// Match <see cref="Key"/> on its <see cref="CvParam.Accession"/> rather than its
    /// <see cref="CvParam.Name"/> — the accession is stable, the display name is not.
    /// PRIDE labels this object <c>"@type": "Tuple"</c> on the wire; that discriminator is ignored,
    /// and the type is named for what it holds rather than for how PRIDE serializes it.
    /// </remarks>
    public class PrideSampleAttribute
    {
        /// <summary>The characteristic being described (e.g. accession "OBI:0100026", name "organism").</summary>
        public CvParam Key { get; set; } = new();

        /// <summary>
        /// The value(s) of that characteristic. Named for the wire field <c>value</c> despite being a
        /// list; not to be confused with <see cref="CvParam.Value"/> on an individual term.
        /// </summary>
        public List<CvParam> Value { get; set; } = new();
    }
}
