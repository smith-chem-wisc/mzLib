using System.Collections.Generic;
using Proteomics;

namespace UsefulProteomicsDatabases.GeneOntology
{
    /// <summary>
    /// Why a protein group's rows say what they say. Declared in ascending precedence: when a group has no
    /// term, the most specific reason for the absence wins, and any term at all makes the group Annotated.
    /// </summary>
    public enum GoAnnotationStatus
    {
        /// <summary>Every member was found in the annotation database, and none carries a GO term.</summary>
        NoGoTerms = 0,

        /// <summary>At least one member is absent from the annotation database (e.g. a FASTA-only accession).</summary>
        NoEntry = 1,

        /// <summary>The group is marked contaminant and has no term. An annotated contaminant reads Annotated.</summary>
        Contaminant = 2,

        /// <summary>
        /// The group has at least one term. A row with a GO id is ALWAYS Annotated -- the invariant consumers
        /// build against -- which is why an annotated contaminant group reads Annotated, not Contaminant.
        /// </summary>
        Annotated = 3
    }

    /// <summary>
    /// A protein group as the annotator needs it: its members and the marks the search gave it. Built by the
    /// caller -- from a stored results file, or from a search's in-memory groups -- so the annotator depends
    /// on no reader and no search engine.
    /// </summary>
    /// <param name="Name">The group's name as the results file writes it, e.g. "P1|P2|P3".</param>
    /// <param name="MemberAccessions">Every member accession. Order is not meaning: MetaMorpheus sorts members
    /// alphabetically and never picks a leading protein.</param>
    /// <param name="IsDecoy">Decoy groups are rejected: decoys carry no GO, so one would pass for a real absence.</param>
    /// <param name="IsContaminant">The search's contaminant mark for the group.</param>
    /// <param name="QValue">The group-level q-value from the results file, carried to every row. Nothing is
    /// filtered on it here: the consumer filters.</param>
    public sealed record GoAnnotationGroup(string Name, IReadOnlyList<string> MemberAccessions, bool IsDecoy,
        bool IsContaminant, double QValue);

    /// <summary>
    /// One row of GO annotation output: a (protein group, GO term) pair, or a group's single term-less row.
    /// Long format: one term per row, never a joined cell of terms. Categories are not here -- they are a
    /// separate per-map table keyed on (go_id, go_release).
    /// </summary>
    /// <param name="ProteinGroup">The group's name.</param>
    /// <param name="AccessionUsed">The members that carry the term, directly or by propagation, in ordinal
    /// order. Empty on a term-less row.</param>
    /// <param name="GoId">The primary GO id, or null on a term-less row.</param>
    /// <param name="GoName">The term's name in the pinned ontology release, or null on a term-less row.</param>
    /// <param name="Aspect">The term's aspect in the pinned release, or null on a term-less row.</param>
    /// <param name="Evidence">ECO codes, in ordinal order, POOLED over every carrying member and every
    /// annotation that produced the term for it: a member's own annotation to the term and its annotations to
    /// descendants, so an evidence filter still bites after propagation. Pooled, so it cannot be combined
    /// with n_with: one member's IDA does not make the others experimental. Use
    /// <paramref name="EvidenceByMember"/> for a per-member filter.</param>
    /// <param name="Inherited">A group-level summary: true when every carrying member is a UniProt isoform (P04406-2) absent from the
    /// annotation database that took the term from its entry (P04406); null on a term-less row. An entry is not
    /// necessarily the isoform's sequence, so the consumer decides whether inherited rows count. Inherited
    /// cellular component terms are usually right, but a location is a positive claim and a wrong one costs
    /// the most: a consumer making organelle claims can rest them on rows where this is false and report the
    /// figure with inherited rows beside it. Nothing here is blanked, and waiting will not fix it: GOA maps
    /// isoform-level subcellular annotation up to the main entry. Which members inherited it is
    /// <paramref name="AccessionInherited"/>.</param>
    /// <param name="Propagated">A group-level summary: true when no carrying member is annotated to the term
    /// itself, only to a descendant; null on a term-less row. Which members are annotated to the term itself
    /// is <paramref name="AccessionDirect"/>. A boolean, never a distance: GO is a DAG, so there is no one path
    /// to measure.</param>
    /// <param name="NMembers">Number of members in the group.</param>
    /// <param name="NWith">Number of members carrying the term; always AccessionUsed.Count.</param>
    /// <param name="Status">The group's annotation status, repeated on each of its rows.</param>
    /// <param name="QValue">The group's q-value, as given.</param>
    /// <param name="GoRelease">The ontology's data-version.</param>
    /// <param name="GoOboSha256">sha256 of the go.obo read.</param>
    /// <param name="AnnotationDbSha256">sha256 of the database the GO terms were read from, as the caller
    /// gave it. It need not be the database that was searched.</param>
    /// <param name="AccessionDirect">The members of <paramref name="AccessionUsed"/> annotated to the term
    /// itself rather than only to a descendant, in ordinal order. Empty on a term-less row. A consensus of
    /// direct annotations is AccessionDirect.Count == NMembers.</param>
    /// <param name="AccessionInherited">The members of <paramref name="AccessionUsed"/> that carry the term
    /// only because they are UniProt isoforms that took their entry's terms, in ordinal order. Empty on a
    /// term-less row.</param>
    /// <param name="EvidenceByMember">Each carrying member's own ECO codes for the term: its annotation to the
    /// term and to any descendant that produced it. Keys are exactly <paramref name="AccessionUsed"/>. Empty
    /// on a term-less row.</param>
    public sealed record GoAnnotationRow(
        string ProteinGroup,
        IReadOnlyList<string> AccessionUsed,
        string GoId,
        string GoName,
        GoAspect? Aspect,
        IReadOnlyList<string> Evidence,
        bool? Inherited,
        bool? Propagated,
        int NMembers,
        int NWith,
        GoAnnotationStatus Status,
        double QValue,
        string GoRelease,
        string GoOboSha256,
        string AnnotationDbSha256,
        IReadOnlyList<string> AccessionDirect,
        IReadOnlyList<string> AccessionInherited,
        IReadOnlyDictionary<string, IReadOnlyList<string>> EvidenceByMember);
}
