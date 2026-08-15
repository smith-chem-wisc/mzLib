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
    internal static class SdrfCell
    {
        /// <summary>
        /// The keys the SDRF-Proteomics specification defines. This is the SOURCE OF TRUTH:
        /// <see cref="KnownLeadingKeys"/> is a subset of <see cref="KnownKeys"/>, pinned by
        /// LeadingKeysAreASubsetOfKnownKeys in the test project, so a key cannot enter the grammar
        /// through the leading-key set without first being justified here.
        ///
        /// NT name, AC accession, MT modification type, TA target amino acid, PP position,
        /// VV version value, CT compound type, QY quantity, PS peptide sequence, SP species,
        /// CN common name, CV vendor, CL cleavable, MH/ML stub mass high/low, SN source name.
        ///
        /// SN is genuinely specified, which is easy to doubt because it is absent from the
        /// modification-style key list: it belongs to characteristics[pooled sample], where
        /// "SN=sample1;SN=sample2" lists the source names of a pool (specification README.adoc
        /// lines 352 and 362, TERMS.tsv "pooled sample"). The corpus agrees -- all 383 leading SN
        /// cells sit in characteristics[pooled sample].
        /// </summary>
        private static readonly HashSet<string> SpecKeys = new(StringComparer.OrdinalIgnoreCase)
        {
            "NT", "AC", "MT", "TA", "PP", "VV", "CT", "QY",
            "PS", "SP", "CN", "CV", "CL", "MH", "ML", "SN"
        };

        /// <summary>
        /// Keys the corpus uses that the specification does not define, kept because they are
        /// widespread and unambiguous in context -- MM monoisotopic mass (26,828 occurrences),
        /// CS cleavage site (3,810), CF chemical formula (2,602), all inside modification and
        /// cleavage-agent cells.
        ///
        /// Separated from <see cref="SpecKeys"/> rather than merged so that the distinction between
        /// "the format says so" and "the community writes it anyway" stays visible. Anything added
        /// here is a tolerance, and tolerances belong on the reading side only -- never in authored
        /// output.
        /// </summary>
        private static readonly HashSet<string> ObservedNonSpecKeys = new(StringComparer.OrdinalIgnoreCase)
        {
            "MM", "CS", "CF"
        };

        /// <summary>
        /// Every key this reader will interpret: the specified ones plus the tolerated ones.
        ///
        /// "PMID" was removed. It is not a key in any position -- the specification uses PMID only
        /// as a VALUE format for reference columns ("URL, DOI, or PMID"), and the 3 corpus
        /// occurrences are of that kind. "N" was removed for the reasons given on
        /// <see cref="KnownLeadingKeys"/>. "PS" was added: it is specified (spiked-compound cells,
        /// "CT=peptide;PS=PEPTIDESEQ;QY=10 fmol") though no corpus file uses it yet.
        /// </summary>
        internal static readonly HashSet<string> KnownKeys =
            new(SpecKeys.Concat(ObservedNonSpecKeys), StringComparer.OrdinalIgnoreCase);

        /// <summary>
        /// Keys that actually appear FIRST in a corpus cell: NT (1,459,508 cells in all 1,236
        /// files), AC (190,369 in 323), CT (389 in 6) and SN (383 in 4). Everything else in
        /// <see cref="KnownKeys"/> only ever appears after a semicolon.
        ///
        /// "A" was removed from the set entirely. It is a single character, it never leads a real
        /// cell, and as a leading key it can only ever produce false positives -- any free-text
        /// value beginning "A=" would have been decoded as a controlled-vocabulary term.
        ///
        /// "N" was removed for exactly the same reasons, and the specification settles it: the key
        /// grammar has no "N", and NT is the only name key. It leads 45 corpus cells -- 27 in
        /// PXD039582 and 18 in PXD039585, every one of them "N=Orbitrap" in
        /// comment[ms2 analyzer type] -- which is a typo for NT= in two files by one submitter, not
        /// a spelling of the grammar. Accepting it decoded drift AS the grammar and so hid it from
        /// the one component whose job is to report it; SdrfDriftLint now sees those cells as the
        /// free text they are. An earlier revision of this comment listed the leading keys as
        /// "NT, AC, CT, SN, and N" with a count against every key except N, which is what made the
        /// reviewer look.
        /// </summary>
        internal static readonly HashSet<string> KnownLeadingKeys = new(StringComparer.OrdinalIgnoreCase)
        {
            "NT", "AC", "CT", "SN", "MT", "TA", "PP", "CS", "CF", "CL", "CV", "MM", "VV"
        };

        /// <summary>
        /// Formats a controlled-vocabulary term as an SDRF cell: "NT=Oxidation;AC=UNIMOD:35".
        ///
        /// NT first, then AC, then any extra keys in the order given. This is a specification MUST,
        /// not merely a convention read off the corpus: "The key order MUST be NT (name) first,
        /// followed by AC (accession), then any additional keys" (README.adoc:261). The corpus
        /// agrees overwhelmingly anyway -- NT leads 1,459,508 cells, AC 190,369.
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

            // NT is the ONLY name key. An earlier revision also accepted "N" as an alternative
            // spelling, on the strength of 45 corpus cells that write "N=Orbitrap"; the
            // specification defines no such key, so that was decoding one submitter's typo as
            // grammar. See KnownLeadingKeys.
            pairs.TryGetValue("NT", out string? name);
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
