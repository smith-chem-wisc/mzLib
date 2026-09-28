using System.Text.RegularExpressions;

namespace MzLibUtil
{
    /// <summary>Which accession grammar an identifier matched.</summary>
    public enum AccessionNamespace
    {
        /// <summary>Matched neither grammar. Kept verbatim; never repaired.</summary>
        Unrecognized = 0,
        UniProt,
        RefSeq
    }

    /// <summary>
    /// An accession as given, plus what can be said about it without any reference data.
    ///
    /// <see cref="EntryAccession"/> is the ENTRY, not "the canonical sequence": P12345-1 is not
    /// necessarily the sequence UniProt displays for P12345, because the displayed isoform can carry
    /// any number. It answers "same UniProt entry?" and nothing about sequence identity. Whether an
    /// isoform counts as the same protein is the consumer's call per question.
    /// </summary>
    /// <param name="Verbatim">Exactly the string given (empty for null).</param>
    /// <param name="EntryAccession">The UniProt entry or unversioned RefSeq accession; Verbatim when unrecognized.</param>
    /// <param name="Isoform">The UniProt isoform number, or null when there is no suffix.</param>
    /// <param name="Version">The RefSeq version, or null when there is none.</param>
    /// <param name="Namespace">Which grammar matched.</param>
    public sealed record ProteinAccession(string Verbatim, string EntryAccession, int? Isoform, int? Version,
        AccessionNamespace Namespace)
    {
        /// <summary>UniProt's published accession grammar, with an optional isoform suffix.</summary>
        private static readonly Regex UniProt = new(
            @"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-(\d+))?$",
            RegexOptions.Compiled);

        /// <summary>RefSeq protein accessions, with an optional version.</summary>
        private static readonly Regex RefSeqProtein = new(@"^((?:NP|XP|YP|WP|AP)_\d+)(?:\.(\d+))?$", RegexOptions.Compiled);

        /// <summary>
        /// Parses, never repairs. Anything outside both grammars -- a decoy or contaminant prefix,
        /// lower case, surrounding whitespace, a FASTA-style "sp|...|..." header -- is Unrecognized
        /// and kept verbatim, so it surfaces as a data problem instead of being silently fixed.
        /// </summary>
        public static ProteinAccession Parse(string accession)
        {
            accession ??= "";

            var m = UniProt.Match(accession);
            if (m.Success)
            {
                return new ProteinAccession(accession, m.Groups[1].Value, NumberOrNull(m.Groups[2]), null,
                    AccessionNamespace.UniProt);
            }

            m = RefSeqProtein.Match(accession);
            if (m.Success)
            {
                return new ProteinAccession(accession, m.Groups[1].Value, null, NumberOrNull(m.Groups[2]),
                    AccessionNamespace.RefSeq);
            }

            return new ProteinAccession(accession, accession, null, null, AccessionNamespace.Unrecognized);
        }

        private static int? NumberOrNull(Group group) =>
            group.Success && int.TryParse(group.Value, out int n) ? n : null;
    }
}
