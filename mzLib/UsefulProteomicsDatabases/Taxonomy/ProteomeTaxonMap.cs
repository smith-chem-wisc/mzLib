#nullable enable
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases.Taxonomy
{
    /// <summary>One UniProt proteome and the taxon it is the proteome of.</summary>
    /// <param name="ProteomeId">The UniProt proteome identifier, e.g. "UP000005640".</param>
    /// <param name="TaxonId">The NCBI taxonomy identifier of the proteome's organism.</param>
    /// <param name="Organism">The organism's name as the source writes it.</param>
    /// <param name="ProteinCount">The number of proteins in the proteome, or null when the source does not say.</param>
    /// <param name="ProteomeType">
    /// The source's own description of the proteome's standing, e.g. "Reference and representative proteome", kept
    /// verbatim; null when the source does not say.
    /// </param>
    public sealed record ProteomeTaxon(
        string ProteomeId,
        long TaxonId,
        string Organism,
        long? ProteinCount,
        string? ProteomeType);

    /// <summary>
    /// Which taxon each UniProt proteome belongs to. A search database downloaded from UniProt is usually named for
    /// its proteome, so this is the step from "the database called UP000005640" to a place in a
    /// <see cref="TaxonomyTree"/>.
    /// </summary>
    public sealed class ProteomeTaxonMap
    {
        internal const string PinnedProteomesResource = "UsefulProteomicsDatabases.Resources.proteome_taxa.tsv.gz";

        private static readonly Lazy<ProteomeTaxonMap> Pinned = new(() =>
            TaxonomySnapshotTsv.ReadProteomesResource(PinnedProteomesResource));

        private readonly Dictionary<string, ProteomeTaxon> _byProteome;
        private readonly List<ProteomeTaxon> _ordered;
        private readonly ILookup<long, ProteomeTaxon> _byTaxon;

        private ProteomeTaxonMap(Dictionary<string, ProteomeTaxon> byProteome, TaxonomyProvenance provenance)
        {
            _byProteome = byProteome;
            Provenance = provenance;
            _ordered = byProteome.Values.OrderBy(p => p.ProteomeId, StringComparer.Ordinal).ToList();
            _byTaxon = _ordered.ToLookup(p => p.TaxonId);
        }

        /// <summary>Where the proteomes came from.</summary>
        public TaxonomyProvenance Provenance { get; }

        /// <summary>How many proteomes the map holds.</summary>
        public int Count => _byProteome.Count;

        /// <summary>Every proteome, in ordinal order of identifier.</summary>
        public IReadOnlyList<ProteomeTaxon> Proteomes => _ordered;

        /// <summary>Builds a map from its proteomes.</summary>
        /// <exception cref="InvalidDataException">A proteome identifier is empty or appears twice.</exception>
        public static ProteomeTaxonMap FromProteomes(IEnumerable<ProteomeTaxon> proteomes, TaxonomyProvenance provenance)
        {
            ArgumentNullException.ThrowIfNull(proteomes);
            ArgumentNullException.ThrowIfNull(provenance);

            var byProteome = new Dictionary<string, ProteomeTaxon>(StringComparer.Ordinal);
            foreach (ProteomeTaxon proteome in proteomes)
            {
                if (string.IsNullOrEmpty(proteome.ProteomeId))
                {
                    throw new InvalidDataException("A proteome has no identifier.");
                }
                if (!byProteome.TryAdd(proteome.ProteomeId, proteome))
                {
                    throw new InvalidDataException($"Proteome {proteome.ProteomeId} appears twice.");
                }
            }
            return new ProteomeTaxonMap(byProteome, provenance);
        }

        /// <summary>
        /// Builds a map from UniProt's downloaded proteome catalogue, as
        /// <see cref="ProteinDbRetriever.UniprotProteomeCatalog"/> reads it. This covers every proteome UniProt
        /// lists, where <see cref="LoadPinned"/> covers the reference proteomes only.
        /// </summary>
        /// <remarks>
        /// A catalogue row with no taxon identifier, or one that is not a whole number, cannot be placed and is left
        /// out; <paramref name="skipped"/> names each one so the caller can report it.
        /// </remarks>
        /// <param name="catalog">The catalogue, keyed by proteome identifier.</param>
        /// <param name="source">Where the catalogue came from, for <see cref="Provenance"/>.</param>
        /// <param name="skipped">The identifiers of the proteomes left out, in ordinal order.</param>
        public static ProteomeTaxonMap FromCatalog(IReadOnlyDictionary<string, UniProtProteomeEntry> catalog, string source,
            out IReadOnlyList<string> skipped)
        {
            ArgumentNullException.ThrowIfNull(catalog);
            ArgumentNullException.ThrowIfNull(source);

            var kept = new List<ProteomeTaxon>();
            var left = new List<string>();
            foreach (UniProtProteomeEntry entry in catalog.Values)
            {
                if (TaxonomyTree.TryParseTaxonId(entry.TaxonId, out long taxonId))
                {
                    kept.Add(new ProteomeTaxon(entry.ProteomeId, taxonId, entry.Organism, entry.ProteinCount, null));
                }
                else
                {
                    left.Add(entry.ProteomeId);
                }
            }
            left.Sort(StringComparer.Ordinal);
            skipped = left;
            return FromProteomes(kept, new TaxonomyProvenance(source));
        }

        /// <summary>
        /// Reads a map from a proteome table written by
        /// <see cref="TaxonomySnapshotTsv.WriteProteomes(TextWriter, ProteomeTaxonMap)"/>. A path ending in ".gz" is
        /// read as gzip.
        /// </summary>
        /// <exception cref="FileNotFoundException">The file does not exist.</exception>
        /// <exception cref="InvalidDataException">The file is not a proteome table.</exception>
        public static ProteomeTaxonMap Load(string proteomesPath) => TaxonomySnapshotTsv.ReadProteomes(proteomesPath);

        /// <summary>
        /// The map shipped inside mzLib: UniProt's reference proteomes, from the same release as
        /// <see cref="TaxonomyTree.LoadPinned"/>, so every taxon it names is in that tree. Read once, on first use.
        /// </summary>
        public static ProteomeTaxonMap LoadPinned() => Pinned.Value;

        /// <summary>Finds a proteome by identifier. The comparison is exact: "UP000005640", not "up000005640".</summary>
        public bool TryGetProteome(string proteomeId, out ProteomeTaxon proteome)
        {
            if (proteomeId == null)
            {
                proteome = null!;
                return false;
            }
            return _byProteome.TryGetValue(proteomeId, out proteome!);
        }

        /// <summary>The taxon a proteome belongs to. False when the map does not hold the proteome.</summary>
        public bool TryGetTaxonId(string proteomeId, out long taxonId)
        {
            bool found = TryGetProteome(proteomeId, out var proteome);
            taxonId = found ? proteome.TaxonId : 0;
            return found;
        }

        /// <summary>
        /// The proteomes of exactly this taxon, in ordinal order of identifier. Empty when there are none. A taxon
        /// can have several: UniProt holds many proteomes for one species.
        /// </summary>
        public IReadOnlyList<ProteomeTaxon> ProteomesOf(long taxonId) => _byTaxon[taxonId].ToList();
    }
}
