using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using MzLibUtil;
using Proteomics;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// Annotates protein groups with Gene Ontology terms without collapsing them.
    ///
    /// For each group it emits one row per GO term held by ANY member -- directly, or by propagation up is_a
    /// and part_of -- and each row records which members carry the term (accession_used) and how many
    /// (n_with of n_members). Which of those members carry it directly, which inherited it, and each one's
    /// evidence are per member too, so a consensus view can be combined with any of them. No member is privileged: MetaMorpheus never chooses a leading protein, and
    /// union, consensus (n_with == n_members) or any other view is a filter the consumer applies to rows.
    ///
    /// Every non-decoy group gets at least one row. A group with no term gets a single term-less row whose
    /// status says why: no member carries GO, a member is absent from the annotation database, or the
    /// group is a contaminant.
    ///
    /// Entrapment is neither target nor decoy, so an entrapment group is annotated, never dropped, and every
    /// row names the group's entrapment members (<see cref="GoAnnotationRow.EntrapmentMembers"/>). A member
    /// is entrapment by the loader's own rule (<see cref="ProteinDbLoader.IsEntrapmentAccession"/>), which
    /// is the only signal a stored MetaMorpheus file carries, or when the annotation database marks its
    /// protein entrapment. An entrapment member contributes its own entry's terms like any other member: a
    /// foreign-proteome entrapment protein's GO is true of its sequence, and the label lets the consumer
    /// decide whether it counts.
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

        /// <summary>Accessions the annotation database marks entrapment, which the accession rule may not see.</summary>
        private readonly HashSet<string> _entrapment = new(StringComparer.Ordinal);

        /// <summary>
        /// The GO ids the annotation database cites that the ontology release does not have, in ordinal order.
        /// Always empty unless the annotator was built with skipUnknownGoIds; pass it to
        /// <see cref="GoAnnotationTsv.Write"/> so the file says terms were dropped.
        /// </summary>
        public IReadOnlyList<string> UnresolvedGoIds { get; }

        /// <param name="ontology">The pinned ontology release terms are resolved and propagated against.</param>
        /// <param name="annotationProteins">Proteins whose GO terms annotate the groups -- typically every target
        /// protein of a UniProt XML. Decoys are ignored. An accession present twice (a target and a contaminant
        /// copy) has its terms unioned.</param>
        /// <param name="annotationDbSha256">sha256 of the database the proteins came from, stamped on every row.</param>
        /// <param name="skipUnknownGoIds">
        /// False (the default) refuses a database that cites a GO id the release lacks. True drops each such id
        /// from its protein and lists it in <see cref="UnresolvedGoIds"/>, so one term newer than a pinned go.obo
        /// does not cost the whole run. A protein whose only terms were dropped then reads no_go_terms.
        /// </param>
        /// <exception cref="ArgumentNullException">A required argument is null.</exception>
        /// <exception cref="InvalidDataException">
        /// Unless <paramref name="skipUnknownGoIds"/>, a protein cites a GO id absent from the ontology release --
        /// usually a UniProt release newer than the pinned go.obo. Every missing id is named. A term that is
        /// merely obsolete is kept.
        /// </exception>
        public GoGroupAnnotator(GeneOntologyGraph ontology, IEnumerable<Protein> annotationProteins, string annotationDbSha256,
            bool skipUnknownGoIds = false)
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
                if (protein.IsEntrapment)
                {
                    _entrapment.Add(protein.Accession);
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

            // term -> (carrying members, members carrying it directly, members carrying it only by inheritance,
            // pooled evidence, and each member's own evidence)
            var byTerm = new SortedDictionary<string, TermAccumulator>(StringComparer.Ordinal);
            var entrapmentMembers = members
                .Where(m => ProteinDbLoader.IsEntrapmentAccession(m) || _entrapment.Contains(m))
                .OrderBy(m => m, StringComparer.Ordinal)
                .ToList();
            bool anyMemberMissing = false;

            foreach (string member in members)
            {
                if (!TryGetTerms(member, out var terms, out bool inherited))
                {
                    anyMemberMissing = true;
                    continue;
                }
                foreach (var (termId, evidence) in terms)
                {
                    Accumulate(byTerm, termId, member, evidence, direct: true, inherited);
                    foreach (string ancestor in _ontology.Ancestors(termId))
                    {
                        Accumulate(byTerm, ancestor, member, evidence, direct: false, inherited);
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
                        _annotationDbSha256, Array.Empty<string>(), Array.Empty<string>(),
                        new Dictionary<string, IReadOnlyList<string>>(StringComparer.Ordinal), entrapmentMembers)
                };
            }

            return byTerm.Select(pair =>
            {
                _ontology.TryGetTerm(pair.Key, out var term);
                var accumulated = pair.Value;
                return new GoAnnotationRow(group.Name, accumulated.Members.ToList(), term.Id, term.Name, term.Aspect,
                    accumulated.Evidence.ToList(),
                    Inherited: accumulated.InheritedMembers.Count == accumulated.Members.Count,
                    Propagated: accumulated.DirectMembers.Count == 0,
                    members.Count, accumulated.Members.Count, GoAnnotationStatus.Annotated, group.QValue,
                    _ontology.Release, _ontology.SourceSha256, _annotationDbSha256,
                    AccessionDirect: accumulated.DirectMembers.ToList(),
                    AccessionInherited: accumulated.InheritedMembers.ToList(),
                    EvidenceByMember: accumulated.EvidenceByMember.ToDictionary(
                        m => m.Key, m => (IReadOnlyList<string>)m.Value.ToList(), StringComparer.Ordinal),
                    EntrapmentMembers: entrapmentMembers);
            }).ToList();
        }

        /// <summary>
        /// Annotates every group, in input order. The rows are materialized, so handing the same result to both
        /// <see cref="GoAnnotationTsv"/> and <see cref="GoCategoryTsv"/> annotates once.
        /// </summary>
        public IReadOnlyList<GoAnnotationRow> AnnotateAll(IEnumerable<GoAnnotationGroup> groups)
        {
            ArgumentNullException.ThrowIfNull(groups);
            return groups.SelectMany(Annotate).ToList();
        }

        /// <summary>
        /// The suffix mzLib appends to a sequence-variant protein's accession
        /// (VariantApplication.GetAccession): one or more "_{original}{position}{variant}" pieces, as in
        /// P04406_A20T or P04406-2_A20T_G31.
        /// </summary>
        private static readonly Regex SequenceVariantSuffix = new(@"^(?<base>[^_]+)(?:_[A-Z*]*\d+[A-Z*]*)+$", RegexOptions.Compiled);

        /// <summary>
        /// A member's own terms when the database has its accession. Otherwise, flagged inherited:
        /// for a sequence variant (P04406_A20T, mzLib's variant accession), the terms of the accession it was
        /// applied to; for a UniProt isoform (P04406-2), its entry's terms; and for a variant of an isoform,
        /// the isoform's or else the entry's. Variants are found when the annotation database was loaded
        /// without variants applied, or is not the one searched. ProteinAccession parses and never repairs, so
        /// an accession outside UniProt's grammar -- a decoy, contaminant or entrapment prefix, a hyphenated
        /// name, a RefSeq NP_ -- never inherits. So an entrapment member absent from the database reads
        /// no_entry rather than borrowing the GO of the target it was made from. Whether an inherited term holds is the consumer's call: an
        /// isoform can differ from its entry precisely in cellular component, which is why the row says so.
        /// </summary>
        private bool TryGetTerms(string member, out Dictionary<string, SortedSet<string>> terms, out bool inherited)
        {
            inherited = false;
            if (_direct.TryGetValue(member, out terms))
            {
                return true;
            }
            inherited = true;
            var variant = SequenceVariantSuffix.Match(member);
            var accession = ProteinAccession.Parse(variant.Success ? variant.Groups["base"].Value : member);
            if (accession.Namespace != AccessionNamespace.UniProt)
            {
                return false;
            }
            if (variant.Success && _direct.TryGetValue(accession.Verbatim, out terms))
            {
                return true;
            }
            return accession.Isoform != null && _direct.TryGetValue(accession.EntryAccession, out terms);
        }

        private static void Accumulate(SortedDictionary<string, TermAccumulator> byTerm, string termId, string member,
            IEnumerable<string> evidence, bool direct, bool inherited)
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
            if (inherited)
            {
                accumulator.InheritedMembers.Add(member);
            }
            accumulator.Evidence.UnionWith(evidence);
            if (!accumulator.EvidenceByMember.TryGetValue(member, out var memberEvidence))
            {
                memberEvidence = new SortedSet<string>(StringComparer.Ordinal);
                accumulator.EvidenceByMember.Add(member, memberEvidence);
            }
            memberEvidence.UnionWith(evidence);
        }

        private sealed class TermAccumulator
        {
            public SortedSet<string> Members { get; } = new(StringComparer.Ordinal);
            public SortedSet<string> DirectMembers { get; } = new(StringComparer.Ordinal);
            public SortedSet<string> InheritedMembers { get; } = new(StringComparer.Ordinal);
            public SortedSet<string> Evidence { get; } = new(StringComparer.Ordinal);
            public SortedDictionary<string, SortedSet<string>> EvidenceByMember { get; } = new(StringComparer.Ordinal);
        }
    }
}
