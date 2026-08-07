using System.Diagnostics.CodeAnalysis;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// Interprets an SDRF cell's key=value grammar, e.g.
    /// "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=Variable".
    ///
    /// This is the opt-in projection that <see cref="SdrfDocument"/> deliberately does not perform
    /// at read time. Cells are stored raw so that documents round-trip losslessly; anything that
    /// wants to reason about a cell's meaning asks for it here.
    ///
    /// Detection is by LEADING KEY, never by shape. A cell is only treated as a term if it starts
    /// with one of the recognised SDRF keys. Testing for "contains '=' and ';'" instead would
    /// misfire on real data: comment[file uri] cells in the curated corpus hold pre-signed download
    /// URLs whose query strings carry "Signature=", "Expires=" and "Id=" -- 5,750 occurrences each
    /// -- and a shape-based parser decodes those as controlled-vocabulary descriptors.
    /// </summary>
    public static class SdrfCell
    {
        /// <summary>
        /// Keys observed across the whole 1,236-file curated corpus, most frequent first:
        /// NT name, AC accession, MT modification type, TA target amino acid, VV version value,
        /// PP position, MM monoisotopic mass, SN specificity name, CS cleavage site,
        /// CF chemical formula, CT compound type, QY quantity, SP species, CN common name,
        /// ML/MH mass low/high, A generic attribute.
        /// </summary>
        private static readonly HashSet<string> KnownKeys = new(StringComparer.OrdinalIgnoreCase)
        {
            "NT", "AC", "MT", "TA", "VV", "PP", "MM", "SN", "CS", "CF",
            "CT", "QY", "SP", "CN", "ML", "MH", "N", "CL", "CV", "PMID"
        };

        /// <summary>
        /// Keys that actually appear FIRST in a corpus cell: NT (1,459,772), AC (190,369),
        /// CT (389), SN (383), and N. Everything else in <see cref="KnownKeys"/> only ever appears
        /// after a semicolon.
        ///
        /// "A" was removed from the set entirely. It is a single character, it never leads a real
        /// cell, and as a leading key it can only ever produce false positives -- any free-text
        /// value beginning "A=" would have been decoded as a controlled-vocabulary term.
        /// </summary>
        private static readonly HashSet<string> KnownLeadingKeys = new(StringComparer.OrdinalIgnoreCase)
        {
            "NT", "AC", "CT", "SN", "N", "MT", "TA", "PP", "CS", "CF", "CL", "CV", "MM", "VV"
        };

        /// <summary>
        /// Formats a controlled-vocabulary term as an SDRF cell: "NT=Oxidation;AC=UNIMOD:35".
        ///
        /// NT first, then AC, then any extra keys in the order given -- the specification's own
        /// examples put the name first, and the curated corpus agrees overwhelmingly (NT leads
        /// 1,459,772 cells, AC 190,369).
        ///
        /// The accession is written EXACTLY as supplied. This deliberately does not "helpfully"
        /// upper-case or add a missing prefix: the corpus is full of drift the caller should not be
        /// able to launder through here -- bare "4" for UNIMOD:4 in 39 documents, "Unimod:35" in 22,
        /// bare "1001251" for Trypsin in 39. Authored terms are meant to come from the pinned
        /// vocabulary already correct; silently repairing a wrong one here would hide the bug and
        /// make the drift lint's job impossible.
        /// </summary>
        /// <param name="term">The term. Name and Accession may not both be empty.</param>
        /// <param name="extras">
        /// Additional key=value pairs in document order, e.g. ("TA","M"), ("MT","Variable").
        /// Keys are emitted as given; a null or empty value is skipped.
        /// </param>
        public static string ToCell(CvParam term, params (string Key, string Value)[] extras)
        {
            if (term is null) throw new ArgumentNullException(nameof(term));
            if (string.IsNullOrEmpty(term.Name) && string.IsNullOrEmpty(term.Accession))
                throw new ArgumentException(
                    "A controlled-vocabulary term needs at least a name or an accession.", nameof(term));

            var parts = new List<string>(2 + (extras?.Length ?? 0));
            if (!string.IsNullOrEmpty(term.Name)) parts.Add("NT=" + term.Name);
            if (!string.IsNullOrEmpty(term.Accession)) parts.Add("AC=" + term.Accession);

            foreach (var (key, value) in extras ?? Array.Empty<(string, string)>())
            {
                if (string.IsNullOrEmpty(key) || string.IsNullOrEmpty(value)) continue;
                parts.Add(key + "=" + value);
            }

            string cell = string.Join(";", parts);

            // A separator inside a value would silently invent or merge keys on read. There is no
            // escape mechanism, so the only honest response is to refuse.
            if (cell.IndexOf('\t') >= 0 || cell.IndexOf('\n') >= 0 || cell.IndexOf('\r') >= 0)
                throw new ArgumentException(
                    $"An SDRF cell cannot contain a tab or newline; the format defines no escape " +
                    $"mechanism. Offending term: '{cell}'.", nameof(term));

            return cell;
        }

        /// <summary>
        /// True when the cell is written in the key=value grammar rather than being free text.
        /// </summary>
        public static bool IsTerm(string cell)
        {
            if (string.IsNullOrEmpty(cell)) return false;
            int equals = cell.IndexOf('=');
            if (equals <= 0) return false;
            int semicolon = cell.IndexOf(';');
            if (semicolon >= 0 && semicolon < equals) return false;

            // Case-insensitive, matching ParseKeyValues. These disagreed before: the leading key was
            // compared Ordinal while every later key was compared OrdinalIgnoreCase, so
            // "NT=Ox;ac=UNIMOD:35" parsed but "nt=Ox;AC=UNIMOD:35" was silently free text.
            return KnownLeadingKeys.Contains(cell.Substring(0, equals).Trim());
        }

        /// <summary>
        /// Splits a cell into its key=value pairs, preserving order and duplicates-last-wins.
        /// Returns an empty dictionary for free text. Keys are upper-cased for lookup; values keep
        /// their original casing and inner whitespace.
        /// </summary>
        public static IReadOnlyDictionary<string, string> ParseKeyValues(string cell)
        {
            var pairs = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            if (!IsTerm(cell)) return pairs;

            foreach (var part in cell.Split(';'))
            {
                int equals = part.IndexOf('=');
                if (equals <= 0) continue;
                string key = part.Substring(0, equals).Trim();
                if (key.Length == 0) continue;
                pairs[key] = part.Substring(equals + 1).Trim();
            }
            return pairs;
        }

        /// <summary>
        /// Reads the cell as a controlled-vocabulary term. False for free text.
        ///
        /// The CV label is derived from the accession prefix ("MS:1001911" -> "MS"), because SDRF
        /// cells do not carry one separately the way an mzML cvParam does.
        /// </summary>
        public static bool TryParseTerm(string cell, [MaybeNullWhen(false)] out CvParam term)
        {
            term = null;
            var pairs = ParseKeyValues(cell);
            if (pairs.Count == 0) return false;

            // "N" is an alternative spelling of the name key -- PXD039582 writes N=Orbitrap in
            // comment[ms2 analyzer type] 27 times.
            if (!pairs.TryGetValue("NT", out string? name))
                pairs.TryGetValue("N", out name);
            pairs.TryGetValue("AC", out string? accession);

            // A cell in the key=value grammar with neither NT nor AC is NOT a controlled-vocabulary
            // term, and must not be promoted to one.
            //
            // An earlier revision fell back to CN/SN/CT/SP here, to stop such cells being reported
            // as free text. That was worse than the problem. ParseKeyValues is last-wins, and the
            // dominant real case is characteristics[pooled sample], where a cell carries up to 45
            // repeated SN= keys ("SN=OSL.53E;SN=OSL.567;..."): the fallback produced a Name that was
            // one arbitrary member of a list, with an empty accession. Worse, callers branch on this
            // method, so promoted cells left the free-text index without entering the accession
            // index -- 538 cells became invisible to every kind of drift analysis at once.
            //
            // Returning false routes them to the free-text side, where their values are still
            // compared. IsTerm still reports true, so a caller that wants the raw pairs can ask
            // ParseKeyValues for them.
            if (string.IsNullOrEmpty(name) && string.IsNullOrEmpty(accession)) return false;

            accession ??= "";
            name ??= "";
            int colon = accession.IndexOf(':');
            string cvLabel = colon > 0 ? accession.Substring(0, colon) : "";

            term = new CvParam(cvLabel, accession, name, "");
            return true;
        }
    }
}
