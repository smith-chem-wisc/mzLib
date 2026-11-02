using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using Proteomics;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// One row of a database's GO annotation: an entry and one GO term it is annotated to directly, or the
    /// entry's single term-less row.
    /// </summary>
    /// <param name="Accession">The database entry's accession, verbatim.</param>
    /// <param name="GoId">The primary GO id (an alternative id cited by UniProt is resolved), or null on a
    /// term-less row: the entry is in the database and carries no GO.</param>
    /// <param name="Aspect">The term's aspect in the pinned release; null on a term-less row.</param>
    /// <param name="Evidence">The entry's own ECO codes for the term, in ordinal order; empty on a term-less row.</param>
    /// <param name="Inherited">True when the entry is a UniProt isoform with no GO of its own and the row is its
    /// entry's term (UniProt writes GO on the entry, not on its isoforms); null on a term-less row.</param>
    public sealed record GoAccessionTermRow(string Accession, string GoId, GoAspect? Aspect,
        IReadOnlyList<string> Evidence, bool? Inherited);

    /// <summary>
    /// The GO annotation of one database, independent of any search: for every entry, the terms it is
    /// annotated to directly, with their evidence. Ancestors are not here -- they are a property of the
    /// ontology release, not of the database (<see cref="GoOntologyTsv.WriteClosure"/>) -- and protein groups
    /// are not here either: which entry a group member takes is <see cref="Annotation.MemberResolver"/>'s.
    ///
    /// So the annotation of a database searched by many deposits is computed and stored once, and a group's
    /// rows are a join: members to entries (the resolver), entries to direct terms (this), direct terms to
    /// ancestors (the closure). <see cref="GoGroupAnnotator"/> computes exactly that join.
    ///
    /// Every entry has at least one row, so "in the database with no GO" stays distinguishable from "not in
    /// the database". A UniProt isoform in the database with no GO of its own takes its entry's terms, flagged
    /// inherited: UniProt writes GO on the entry, so a database holding isoforms as separate entries holds
    /// every isoform with none, and reading that as "no GO" would make the answer depend on which database
    /// was used. An isoform with GO of its own keeps it.
    /// </summary>
    public sealed class GoAccessionAnnotation
    {
        /// <summary>accession -> primary GO id -> the entry's own evidence codes for it.</summary>
        private readonly Dictionary<string, Dictionary<string, SortedSet<string>>> _direct = new(StringComparer.Ordinal);

        /// <param name="ontology">The pinned ontology release terms are resolved against.</param>
        /// <param name="annotationProteins">The database's proteins, typically every target protein of a UniProt
        /// XML. Decoys are ignored. An accession present twice (a target and a contaminant copy) has its terms
        /// unioned.</param>
        /// <param name="annotationDbSha256">sha256 of the database the proteins came from.</param>
        /// <param name="skipUnknownGoIds">
        /// False (the default) refuses a database that cites a GO id the release lacks. True drops each such id
        /// from its protein and lists it in <see cref="UnresolvedGoIds"/>, so one term newer than a pinned go.obo
        /// does not cost the whole run. A protein whose only terms were dropped then has no GO.
        /// </param>
        /// <exception cref="ArgumentNullException">ontology or annotationProteins is null.</exception>
        /// <exception cref="InvalidDataException">
        /// Unless <paramref name="skipUnknownGoIds"/>, a protein cites a GO id absent from the ontology release --
        /// usually a UniProt release newer than the pinned go.obo. Every missing id is named. A term that is
        /// merely obsolete is kept.
        /// </exception>
        public GoAccessionAnnotation(GeneOntologyGraph ontology, IEnumerable<Protein> annotationProteins,
            string annotationDbSha256, bool skipUnknownGoIds = false)
        {
            Ontology = ontology ?? throw new ArgumentNullException(nameof(ontology));
            ArgumentNullException.ThrowIfNull(annotationProteins);
            AnnotationDbSha256 = annotationDbSha256;

            var missing = new SortedSet<string>(StringComparer.Ordinal);
            foreach (var protein in annotationProteins)
            {
                if (protein == null || protein.IsDecoy)
                {
                    continue;
                }
                if (!_direct.TryGetValue(protein.Accession, out var terms))
                {
                    terms = new Dictionary<string, SortedSet<string>>(StringComparer.Ordinal);
                    _direct.Add(protein.Accession, terms);
                }
                foreach (var goTerm in protein.GoTerms)
                {
                    if (!ontology.TryGetTerm(goTerm.Id, out var term))
                    {
                        missing.Add(goTerm.Id);
                        continue;
                    }
                    if (!terms.TryGetValue(term.Id, out var evidence))
                    {
                        evidence = new SortedSet<string>(StringComparer.Ordinal);
                        terms.Add(term.Id, evidence);
                    }
                    evidence.UnionWith(goTerm.EvidenceCodes);
                }
            }

            if (missing.Count > 0 && !skipUnknownGoIds)
            {
                throw new InvalidDataException(
                    $"The annotation database cites {missing.Count} GO id(s) absent from Gene Ontology release " +
                    $"{ontology.Release ?? "(unversioned)"}: {string.Join(", ", missing)}. Load a go.obo release at least " +
                    "as new as the database (Loaders.UpdateGeneOntology), or skip unknown ids and report them.");
            }
            UnresolvedGoIds = missing.ToList();
            Accessions = _direct.Keys.OrderBy(a => a, StringComparer.Ordinal).ToList();
        }

        /// <summary>The release terms were resolved against.</summary>
        public GeneOntologyGraph Ontology { get; }

        /// <summary>sha256 of the annotation database, as the caller gave it.</summary>
        public string AnnotationDbSha256 { get; }

        /// <summary>
        /// The GO ids the database cites that the release does not have, in ordinal order. Always empty unless
        /// built with skipUnknownGoIds.
        /// </summary>
        public IReadOnlyList<string> UnresolvedGoIds { get; }

        /// <summary>Every entry of the database, in ordinal order.</summary>
        public IReadOnlyList<string> Accessions { get; }

        /// <summary>
        /// An entry's direct terms (primary id -> its evidence codes) and whether they are its entry's, taken by
        /// a UniProt isoform with no GO of its own. False when the database has no such entry.
        /// </summary>
        public bool TryGetTerms(string accession, out IReadOnlyDictionary<string, SortedSet<string>> terms, out bool inherited)
        {
            inherited = false;
            terms = null;
            if (accession == null || !_direct.TryGetValue(accession, out var own))
            {
                return false;
            }
            terms = own;
            if (own.Count == 0
                && ProteinAccession.Parse(accession) is { Namespace: AccessionNamespace.UniProt, Isoform: not null } isoform
                && _direct.TryGetValue(isoform.EntryAccession, out var entryTerms)
                && entryTerms.Count > 0)
            {
                terms = entryTerms;
                inherited = true;
            }
            return true;
        }

        /// <summary>
        /// Every row: entries in ordinal order, terms within an entry in ordinal order, and one term-less row for
        /// an entry with no term.
        /// </summary>
        public IReadOnlyList<GoAccessionTermRow> Rows()
        {
            var rows = new List<GoAccessionTermRow>();
            foreach (string accession in Accessions)
            {
                TryGetTerms(accession, out var terms, out bool inherited);
                if (terms.Count == 0)
                {
                    rows.Add(new GoAccessionTermRow(accession, null, null, Array.Empty<string>(), null));
                    continue;
                }
                foreach (var (goId, evidence) in terms.OrderBy(t => t.Key, StringComparer.Ordinal))
                {
                    Ontology.TryGetTerm(goId, out var term);
                    rows.Add(new GoAccessionTermRow(accession, goId, term.Aspect, evidence.ToList(), inherited));
                }
            }
            return rows;
        }
    }
}
