using System.Collections.Generic;
using MzLibUtil;
using Newtonsoft.Json.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One dataset found by <see cref="MassiveClient.SearchDatasetsAsync"/>: what MassIVE's search shows for it. Use
    /// <see cref="MassiveClient.GetDatasetAsync"/> for the full record.
    /// </summary>
    /// <remarks>
    /// Deliberately not a <see cref="MassiveDataset"/>. A search row leaves fields out (no summary, keywords,
    /// modifications or links), and a <see cref="MassiveDataset"/> with those null would claim the record does not
    /// state them. The fields it has follow the record's null rule: null means the row does not state it.
    /// Contacts are not exposed: they carry email addresses.
    /// </remarks>
    public sealed class MassiveDatasetSearchResult
    {
        /// <summary>The MassIVE accession, e.g. "MSV000078304". Never null.</summary>
        public string MassiveAccession { get; init; }

        /// <summary>The ProteomeXchange accession, or null when the dataset has none (48 of the 67 "zebrafish" hits, 2026-11-10).</summary>
        public string PxdAccession { get; init; }

        /// <summary>The title, with a leading "Title: " cut, or null when not stated.</summary>
        public string Title { get; init; }

        /// <summary>The species, one entry per term group, or null when not stated.</summary>
        public IReadOnlyList<MassiveSpecies> Species { get; init; }

        /// <summary>The instruments, or null when not stated; as <see cref="MassiveDataset.Instruments"/>.</summary>
        public IReadOnlyList<CvParam> Instruments { get; init; }

        /// <summary>The publications, or null when not stated; as <see cref="MassiveDataset.Publications"/>.</summary>
        public IReadOnlyList<MassivePublication> Publications { get; init; }

        /// <summary>One compact PROXI row, or null when it names no MassIVE accession.</summary>
        internal static MassiveDatasetSearchResult FromProxi(JObject row)
        {
            (string msv, string pxd) = MassiveDataset.AccessionsOf(row["accession"]);
            if (msv == null)
                return null;

            return new MassiveDatasetSearchResult
            {
                MassiveAccession = msv,
                PxdAccession = pxd,
                Title = MassiveDataset.TitleOf(row["title"]).Title,
                Species = MassiveDataset.SpeciesOf(row["species"]),
                Instruments = MassiveDataset.InstrumentsOf(row["instruments"]),
                Publications = MassiveDataset.PublicationsOf(row["publications"]),
            };
        }

        /// <inheritdoc />
        public override string ToString() => $"{MassiveAccession}: {Title}";
    }
}
