using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Proteomics;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// Annotates protein groups with Gene Ontology terms without collapsing them.
    ///
    /// For each group it emits one row per GO term held by ANY member -- directly, or by propagation up is_a
    /// and part_of -- and each row records which members carry the term (accession_used) and how many
    /// (n_with of n_members). No member is privileged: MetaMorpheus never chooses a leading protein, and
    /// union, consensus (n_with == n_members) or any other view is a filter the consumer applies to rows.
    ///
    /// Every non-decoy group gets at least one row. A group with no term gets a single term-less row whose
    /// status says why: no member carries GO, a member is absent from the annotation database, or the
    /// group is a contaminant.
    ///
    /// The annotator takes proteins already loaded (ProteinDbLoader.LoadProteinXML mutates static state and
    /// is not safe to call from here) and never filters on q-value or evidence: both travel on every row.
    /// </summary>
    public sealed class GoGroupAnnotator
    {
        private readonly GeneOntologyGraph _ontology;
        private readonly string _annotationDbSha256;

        /// <summary>accession -> primary GO id -> evidence codes of the member's own annotation to that term.</summary>
        private readonly Dictionary<string, Dictionary<string, SortedSet<string>>> _direct = new(StringComparer.Ordinal);

        /// <param name="ontology">The pinned ontology release terms are resolved and propagated against.</param>
        /// <param name="annotationProteins">Proteins whose GO terms annotate the groups -- typically every target
        /// protein of a UniProt XML. Decoys are ignored. An accession present twice (a target and a contaminant
        /// copy) has its terms unioned.</param>
        /// <param name="annotationDbSha256">sha256 of the database the proteins came from, stamped on every row.</param>
        /// <exception cref="ArgumentNullException">A required argument is null.</exception>
        /// <exception cref="InvalidDataException">
        /// A protein cites a GO id absent from the ontology release -- usually a UniProt release newer than the
        /// pinned go.obo. Every missing id is named. A term that is merely obsolete is kept.
        /// </exception>
        public GoGroupAnnotator(GeneOntologyGraph ontology, IEnumerable<Protein> annotationProteins, string annotationDbSha256)
        {
            _ontology = ontology ?? throw new ArgumentNullException(nameof(ontology));
            ArgumentNullException.ThrowIfNull(annotationProteins);
            _annotationDbSha256 = annotationDbSha256;

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

            if (missing.Count > 0)
            {
                throw new InvalidDataException(
                    $"The annotation database cites {missing.Count} GO id(s) absent from Gene Ontology release " +
                    $"{ontology.Release ?? "(unversioned)"}: {string.Join(", ", missing)}. Load a go.obo release at least as new as the database.");
            }
        }

        /// <summary>
        /// The rows for one group, ordered by GO id (ordinal); a term-less group gets exactly one row.
        /// </summary>
        /// <exception cref="ArgumentException">The group is a decoy, or has no members.</exception>
        public IReadOnlyList<GoAnnotationRow> Annotate(GoAnnotationGroup group)
        {
            ArgumentNullException.ThrowIfNull(group);
            if (group.IsDecoy)
            {
                throw new ArgumentException($"Decoy group '{group.Name}' cannot be annotated: decoys carry no GO.", nameof(group));
            }
            var members = (group.MemberAccessions ?? Array.Empty<string>())
                .Where(a => !string.IsNullOrEmpty(a))
                .Distinct(StringComparer.Ordinal)
                .ToList();
            if (members.Count == 0)
            {
                throw new ArgumentException($"Group '{group.Name}' has no members.", nameof(group));
            }

            // term -> (carrying members, members carrying it directly, evidence)
            var byTerm = new SortedDictionary<string, TermAccumulator>(StringComparer.Ordinal);
            bool anyMemberMissing = false;

            foreach (string member in members)
            {
                if (!_direct.TryGetValue(member, out var terms))
                {
                    anyMemberMissing = true;
                    continue;
                }
                foreach (var (termId, evidence) in terms)
                {
                    Accumulate(byTerm, termId, member, evidence, direct: true);
                    foreach (string ancestor in _ontology.Ancestors(termId))
                    {
                        Accumulate(byTerm, ancestor, member, evidence, direct: false);
                    }
                }
            }

            if (byTerm.Count == 0)
            {
                var status = group.IsContaminant ? GoAnnotationStatus.Contaminant
                    : anyMemberMissing ? GoAnnotationStatus.NoEntry
                    : GoAnnotationStatus.NoGoTerms;
                return new[]
                {
                    new GoAnnotationRow(group.Name, Array.Empty<string>(), null, null, null, Array.Empty<string>(),
                        null, null, members.Count, 0, status, group.QValue, _ontology.Release, _ontology.SourceSha256,
                        _annotationDbSha256)
                };
            }

            return byTerm.Select(pair =>
            {
                _ontology.TryGetTerm(pair.Key, out var term);
                var accumulated = pair.Value;
                return new GoAnnotationRow(group.Name, accumulated.Members.ToList(), term.Id, term.Name, term.Aspect,
                    accumulated.Evidence.ToList(),
                    Inherited: false,
                    Propagated: accumulated.DirectMembers.Count == 0,
                    members.Count, accumulated.Members.Count, GoAnnotationStatus.Annotated, group.QValue,
                    _ontology.Release, _ontology.SourceSha256, _annotationDbSha256);
            }).ToList();
        }

        /// <summary>Annotates every group, in input order.</summary>
        public IEnumerable<GoAnnotationRow> AnnotateAll(IEnumerable<GoAnnotationGroup> groups)
        {
            ArgumentNullException.ThrowIfNull(groups);
            return groups.SelectMany(Annotate);
        }

        private static void Accumulate(SortedDictionary<string, TermAccumulator> byTerm, string termId, string member,
            IEnumerable<string> evidence, bool direct)
        {
            if (!byTerm.TryGetValue(termId, out var accumulator))
            {
                accumulator = new TermAccumulator();
                byTerm.Add(termId, accumulator);
            }
            accumulator.Members.Add(member);
            if (direct)
            {
                accumulator.DirectMembers.Add(member);
            }
            accumulator.Evidence.UnionWith(evidence);
        }

        private sealed class TermAccumulator
        {
            public SortedSet<string> Members { get; } = new(StringComparer.Ordinal);
            public SortedSet<string> DirectMembers { get; } = new(StringComparer.Ordinal);
            public SortedSet<string> Evidence { get; } = new(StringComparer.Ordinal);
        }
    }
}
