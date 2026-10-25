using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using MzLibUtil;
using Proteomics;

namespace UsefulProteomicsDatabases.Annotation
{
    /// <summary>How a protein-group member was matched to an entry of the annotation database.</summary>
    public enum MemberResolutionKind
    {
        /// <summary>The database has the member's own accession.</summary>
        Self = 0,

        /// <summary>
        /// The database lacks the member, and the member takes another entry's annotation by UniProt's or mzLib's
        /// accession grammar: a UniProt isoform (P04406-2) its entry's (P04406), a sequence variant (P04406_A20T)
        /// the accession it was applied to.
        /// </summary>
        Inherited = 1,

        /// <summary>The database lacks the member and the grammar gives it no other entry.</summary>
        NoEntry = 2
    }

    /// <summary>
    /// One member accession of a search's protein groups and the database entry whose annotation it takes.
    /// </summary>
    /// <param name="Accession">The member exactly as the search wrote it.</param>
    /// <param name="AnnotatedAs">The database accession whose annotation rows the member takes: itself, the entry
    /// of a UniProt isoform, or the accession a variant was applied to. Null when <paramref name="Kind"/> is
    /// <see cref="MemberResolutionKind.NoEntry"/>.</param>
    /// <param name="Kind">How the member was matched.</param>
    /// <param name="IsEntrapment">The member is entrapment, by the loader's accession rule or because the
    /// database marks it. Entrapment is neither target nor decoy: it is resolved like any member and labelled.</param>
    /// <param name="SequenceAccession">The sequence the member is a variant of, or the member itself
    /// (<see cref="MemberResolver.SequenceAccessionOf"/>): P04406-2_A20T is isoform 2 with a substitution, so
    /// P04406-2, whichever entry it takes its annotation from. An annotation that names one isoform applies to the
    /// members whose sequence accession is that isoform. Null when <paramref name="Kind"/> is
    /// <see cref="MemberResolutionKind.NoEntry"/>.</param>
    public sealed record MemberResolution(string Accession, string AnnotatedAs, MemberResolutionKind Kind, bool IsEntrapment,
        string SequenceAccession);

    /// <summary>
    /// Matches the members of a search's protein groups to the entries of an annotation database.
    ///
    /// A group's members are not all database entries. A search can report a UniProt isoform the database
    /// holds only as its entry, or a sequence variant it built at search time (mzLib writes P04406_A20T). This
    /// says, for each member, which database accession's annotation it takes. It depends only on which
    /// accessions the database holds, never on what they are annotated with, so one resolution serves any
    /// annotation read from that database: GO, pathways, or anything else carried by a dbReference.
    ///
    /// The rule parses and never repairs (<see cref="ProteinAccession.Parse"/>): an accession outside UniProt's
    /// grammar -- a decoy, contaminant or entrapment prefix, a hyphenated name, a RefSeq NP_ -- never takes
    /// another entry's annotation. So an entrapment member absent from the database reads NoEntry rather
    /// than borrowing the annotation of the target it was made from.
    /// </summary>
    public sealed class MemberResolver
    {
        /// <summary>
        /// The suffix mzLib appends to a sequence-variant protein's accession
        /// (VariantApplication.GetAccession): one or more "_{original}{position}{variant}" pieces, as in
        /// P04406_A20T or P04406-2_A20T_G31.
        /// </summary>
        private static readonly Regex SequenceVariantSuffix = new(@"^(?<base>[^_]+)(?:_[A-Z*]*\d+[A-Z*]*)+$", RegexOptions.Compiled);

        private readonly HashSet<string> _accessions = new(StringComparer.Ordinal);
        private readonly HashSet<string> _entrapment = new(StringComparer.Ordinal);

        /// <param name="databaseProteins">The annotation database's proteins, typically every target protein of
        /// the searched UniProt XML. Decoys are ignored: they are not entries.</param>
        /// <exception cref="ArgumentNullException">databaseProteins is null.</exception>
        public MemberResolver(IEnumerable<Protein> databaseProteins)
        {
            ArgumentNullException.ThrowIfNull(databaseProteins);
            foreach (var protein in databaseProteins)
            {
                if (protein == null || protein.IsDecoy)
                {
                    continue;
                }
                _accessions.Add(protein.Accession);
                if (protein.IsEntrapment)
                {
                    _entrapment.Add(protein.Accession);
                }
            }
        }

        /// <summary>True when the database holds this accession as an entry.</summary>
        public bool Contains(string accession) => accession != null && _accessions.Contains(accession);

        /// <summary>True when the member is entrapment by the loader's accession rule or the database's mark.</summary>
        public bool IsEntrapment(string accession) =>
            accession != null && (ProteinDbLoader.IsEntrapmentAccession(accession) || _entrapment.Contains(accession));

        /// <summary>
        /// Resolves one member. A member in the database is itself. Otherwise a sequence variant takes the
        /// accession it was applied to when the database has it, and a UniProt isoform -- or a variant of one --
        /// takes its entry.
        /// </summary>
        /// <exception cref="ArgumentException">The accession is null or empty.</exception>
        public MemberResolution Resolve(string accession)
        {
            if (string.IsNullOrEmpty(accession))
            {
                throw new ArgumentException("A member accession is null or empty.", nameof(accession));
            }
            bool entrapment = IsEntrapment(accession);
            var variant = SequenceVariantSuffix.Match(accession);
            var parsed = ProteinAccession.Parse(variant.Success ? variant.Groups["base"].Value : accession);
            bool uniProt = parsed.Namespace == AccessionNamespace.UniProt;
            string sequence = variant.Success && uniProt ? parsed.Verbatim : accession;

            if (_accessions.Contains(accession))
            {
                return new MemberResolution(accession, accession, MemberResolutionKind.Self, entrapment, sequence);
            }
            if (uniProt)
            {
                if (variant.Success && _accessions.Contains(parsed.Verbatim))
                {
                    return new MemberResolution(accession, parsed.Verbatim, MemberResolutionKind.Inherited, entrapment, sequence);
                }
                if (parsed.Isoform != null && _accessions.Contains(parsed.EntryAccession))
                {
                    return new MemberResolution(accession, parsed.EntryAccession, MemberResolutionKind.Inherited, entrapment, sequence);
                }
            }
            return new MemberResolution(accession, null, MemberResolutionKind.NoEntry, entrapment, null);
        }

        /// <summary>
        /// The sequence an accession names: with mzLib's sequence-variant suffix removed when what remains is a
        /// UniProt accession (P04406-2_A20T is P04406-2, P04406_A20T is P04406), and otherwise the accession itself.
        /// The UniProt condition keeps an accession that merely looks like a variant (CON_P02768) whole.
        /// </summary>
        /// <exception cref="ArgumentException">The accession is null or empty.</exception>
        public static string SequenceAccessionOf(string accession)
        {
            if (string.IsNullOrEmpty(accession))
            {
                throw new ArgumentException("An accession is null or empty.", nameof(accession));
            }
            var variant = SequenceVariantSuffix.Match(accession);
            if (!variant.Success)
            {
                return accession;
            }
            var parsed = ProteinAccession.Parse(variant.Groups["base"].Value);
            return parsed.Namespace == AccessionNamespace.UniProt ? parsed.Verbatim : accession;
        }

        /// <summary>
        /// Resolves every distinct member, in ordinal order: the per-deposit table, one row per accession. Null
        /// and empty accessions are skipped.
        /// </summary>
        /// <exception cref="ArgumentNullException">accessions is null.</exception>
        public IReadOnlyList<MemberResolution> ResolveAll(IEnumerable<string> accessions)
        {
            ArgumentNullException.ThrowIfNull(accessions);
            return accessions
                .Where(a => !string.IsNullOrEmpty(a))
                .Distinct(StringComparer.Ordinal)
                .OrderBy(a => a, StringComparer.Ordinal)
                .Select(Resolve)
                .ToList();
        }
    }
}
