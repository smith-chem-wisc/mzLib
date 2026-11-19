using System;
using System.Collections.Generic;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// One list from ProteomeCentral's search and MassIVE's own: each dataset once, keyed by its PXD, else its MSV, and
    /// each saying which search found it (aging, PRIDE-M12).
    /// </summary>
    /// <remarks>
    /// The two searches match different fields, so a dataset's presence means different things. ProteomeCentral's
    /// <c>search=</c> matches titles and submitter keywords; MassIVE's matches inside words ("imaging" for "aging") and
    /// reaches descriptions (aging's survey, 2026-11-10). And only MassIVE's sees the datasets that have no PXD. So each
    /// hit keeps both sources' rows as they came, rather than one merged row that would hide which said what.
    /// </remarks>
    public static class RepositorySearch
    {
        /// <summary>
        /// Merges a ProteomeCentral search with a MassIVE one: every ProteomeCentral hit in its order, joined by PXD with
        /// the MassIVE hit for the same dataset, then the MassIVE hits ProteomeCentral did not return, in MassIVE's order.
        /// </summary>
        /// <param name="central">ProteomeCentral's hits, or null for none.</param>
        /// <param name="massive">MassIVE's hits, or null for none.</param>
        /// <returns>Each dataset once. Never null.</returns>
        public static IReadOnlyList<RepositorySearchHit> Merge(IEnumerable<ProteomeCentralDataset> central,
            IEnumerable<MassiveDatasetSearchResult> massive)
        {
            var hits = new List<RepositorySearchHit>();
            var byPxd = new Dictionary<string, int>(StringComparer.OrdinalIgnoreCase);

            foreach (ProteomeCentralDataset row in central ?? Array.Empty<ProteomeCentralDataset>())
            {
                if (row?.Accession == null || byPxd.ContainsKey(row.Accession))
                    continue;
                byPxd[row.Accession] = hits.Count;
                hits.Add(new RepositorySearchHit(row, null));
            }

            var byMsv = new HashSet<string>(StringComparer.OrdinalIgnoreCase);
            foreach (MassiveDatasetSearchResult row in massive ?? Array.Empty<MassiveDatasetSearchResult>())
            {
                if (row?.MassiveAccession == null || !byMsv.Add(row.MassiveAccession))
                    continue;
                if (row.PxdAccession != null && byPxd.TryGetValue(row.PxdAccession, out int at))
                {
                    // The first MassIVE row for a PXD is joined; a second one (another MSV under the same PXD) adds nothing.
                    if (hits[at].Massive == null)
                        hits[at] = new RepositorySearchHit(hits[at].Central, row);
                }
                else
                {
                    if (row.PxdAccession != null)
                        byPxd[row.PxdAccession] = hits.Count;
                    hits.Add(new RepositorySearchHit(null, row));
                }
            }

            return hits;
        }
    }

    /// <summary>One dataset of a <see cref="RepositorySearch.Merge"/>: what each search said about it.</summary>
    public sealed class RepositorySearchHit
    {
        internal RepositorySearchHit(ProteomeCentralDataset central, MassiveDatasetSearchResult massive)
        {
            Central = central;
            Massive = massive;
        }

        /// <summary>The PXD, from either search; null for a MassIVE dataset that has none.</summary>
        public string PxdAccession => Central?.Accession ?? Massive?.PxdAccession;

        /// <summary>The MSV, when MassIVE's search returned the dataset (ProteomeCentral's rows do not carry it).</summary>
        public string MassiveAccession => Massive?.MassiveAccession;

        /// <summary>The title: ProteomeCentral's when it found the dataset, else MassIVE's.</summary>
        public string Title => Central?.Title ?? Massive?.Title;

        /// <summary>Whether ProteomeCentral's search returned the dataset.</summary>
        public bool FoundByProteomeCentral => Central != null;

        /// <summary>Whether MassIVE's search returned the dataset.</summary>
        public bool FoundByMassive => Massive != null;

        /// <summary>ProteomeCentral's row, as it came; null when its search did not return the dataset.</summary>
        public ProteomeCentralDataset Central { get; }

        /// <summary>MassIVE's row, as it came; null when its search did not return the dataset.</summary>
        public MassiveDatasetSearchResult Massive { get; }

        /// <inheritdoc />
        public override string ToString() =>
            $"{PxdAccession ?? MassiveAccession} [{(FoundByProteomeCentral ? "proteomecentral" : "")}{(FoundByProteomeCentral && FoundByMassive ? "+" : "")}{(FoundByMassive ? "massive" : "")}]: {Title}";
    }
}
