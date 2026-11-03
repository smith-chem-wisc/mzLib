using System;
using System.Collections.Generic;
using Proteomics;

namespace UsefulProteomicsDatabases.Pathways
{
    /// <summary>What one entry's Reactome rows say about it.</summary>
    public enum ReactomeMembershipStatus
    {
        /// <summary>The row carries a pathway, the entry's own or inherited.</summary>
        Annotated = 0,

        /// <summary>
        /// UniProt lists no Reactome reference for the entry, and it inherits none, while its organism is covered.
        /// This means "UniProt lists none", not "in no pathway".
        /// </summary>
        NoReference = 1,

        /// <summary>No entry of the entry's organism in the loaded list carries any Reactome reference.</summary>
        SpeciesNotCovered = 2
    }

    /// <summary>
    /// One row of a database's Reactome membership: an entry and one pathway it is a member of, or the entry's single
    /// pathway-less row.
    /// </summary>
    /// <param name="Accession">The database entry's accession, verbatim.</param>
    /// <param name="Organism">The entry's own NCBI taxon as a CURIE, "NCBITaxon:9606"; null when the entry names none.</param>
    /// <param name="PathwayId">The Reactome stable id, or null on a pathway-less row.</param>
    /// <param name="Molecule">The isoform the membership is scoped to (UniProt's &lt;molecule&gt;), or null when it covers
    /// the whole entry.</param>
    /// <param name="SequenceAccession">Null on a whole-entry row. On an isoform row: the entry's own sequence accession
    /// (<see cref="Annotation.MemberResolver.SequenceAccessionOf"/>) when the molecule is the entry's displayed isoform,
    /// otherwise the molecule. A member takes the isoform rows whose sequence_accession is its own.</param>
    /// <param name="Inherited">True when the entry is a UniProt isoform with no reference of its own and the row is
    /// its entry's whole-entry row; null on a pathway-less row.</param>
    /// <param name="Status">Annotated on every row that carries a pathway.</param>
    public sealed record ReactomeMembershipRow(string Accession, string Organism, string PathwayId, string Molecule,
        string SequenceAccession, bool? Inherited, ReactomeMembershipStatus Status);

    /// <summary>
    /// The Reactome membership of one annotation database, independent of any search
    /// (kegg:DEF-REACTOME-MEMBERSHIP v1). For every entry: the pathways UniProt lists on it, each scoped to the whole
    /// entry or to one isoform, after a species check. Ancestors are not here; they are a property of the Reactome
    /// release (the per-release Reactome tables). Which entry a protein-group member takes is
    /// <see cref="Annotation.MemberResolver"/>'s.
    /// </summary>
    public sealed class ReactomeAccessionAnnotation
    {
        /// <summary>The definition every row follows, written in the header.</summary>
        public const string DefinitionId = "kegg:DEF-REACTOME-MEMBERSHIP";

        /// <summary>The definition's version. It changes only when the same inputs would give different rows.</summary>
        public const int DefinitionVersion = 1;

        /// <summary>An organism needs at least this many Reactome references for its species code to be read.</summary>
        public const int MinimumReferencesForSpeciesCode = 10;

        /// <summary>
        /// The name the species rule and the header give the entries that name no NCBI taxon (a FASTA entry without
        /// OX=). They are treated as one organism; their rows' organism is null.
        /// </summary>
        public const string UnknownOrganism = "unknown";

        /// <param name="annotationProteins">The target proteins of one search, every file it loaded. Decoys and nulls are
        /// skipped. An accession present more than once has its references unioned.</param>
        /// <param name="annotationDbSha256">The target-set hash the artefact is keyed by.</param>
        /// <exception cref="ArgumentNullException">annotationProteins is null.</exception>
        /// <exception cref="ArgumentException">annotationDbSha256 is null or empty: the table is per database.</exception>
        public ReactomeAccessionAnnotation(IEnumerable<Protein> annotationProteins, string annotationDbSha256)
        {
            throw new NotImplementedException();
        }

        /// <summary>The target-set hash, as the caller gave it.</summary>
        public string AnnotationDbSha256 { get; }

        /// <summary>Every entry, in ordinal order.</summary>
        public IReadOnlyList<string> Accessions { get; } = Array.Empty<string>();

        /// <summary>Each organism whose code was read (a strict majority of at least
        /// <see cref="MinimumReferencesForSpeciesCode"/> references), and its code, e.g. "HSA".</summary>
        public IReadOnlyDictionary<string, string> SpeciesCodes { get; } = new Dictionary<string, string>();

        /// <summary>Each organism with references but no code, and its reference count. Its references are kept
        /// unfiltered.</summary>
        public IReadOnlyDictionary<string, int> UnvalidatedOrganisms { get; } = new Dictionary<string, int>();

        /// <summary>The references dropped because their prefix is not their organism's code.</summary>
        public int DroppedForeignReferences { get; }

        /// <summary>Every row: entries in ordinal order, then pathway id, then molecule (whole-entry first).</summary>
        public IReadOnlyList<ReactomeMembershipRow> Rows() => throw new NotImplementedException();
    }
}
