using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using MzLibUtil;
using TopDownProteomics.IO.Obo;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A pinned, embedded controlled vocabulary — PSI-MS or the PRIDE CV — with lookup by accession
    /// and by name.
    ///
    /// EMBEDDED AND PINNED, deliberately, rather than resolved against a live service. Resolving a
    /// term against something that changes means the same experiment annotated a year apart gets
    /// different accessions as the vocabulary drifts, and a corpus built that way silently stops
    /// joining to itself. The version each vocabulary was built from is exposed as
    /// <see cref="Version"/> so it can be recorded in whatever the terms end up in; updating a
    /// vocabulary is then a deliberate, reviewable change to this repository rather than something
    /// that happens to a user overnight.
    ///
    /// Size is unremarkable next to what mzLib already ships: psi-ms.obo is 1.1 MB and pride_cv.obo
    /// 224 KB, against the 4.4 MB PSI-MOD.obo.xml and 2.3 MB unimod.xml already embedded in Omics.
    /// The large sample ontologies — NCBITaxon at roughly 600 MB, UBERON, MONDO — are deliberately
    /// NOT here: an organism identifier comes from the search database that was already loaded
    /// (see ProteinDbLoader.NcbiTaxonomyDatabaseReferenceType), so it never has to be looked up.
    ///
    /// Parsing reuses TopDownProteomics' OboParser, already a dependency of this project and
    /// already used by <see cref="Loaders.ReadPsiModFile"/>. Loading is lazy, matching
    /// ProteaseDictionary's cached-embedded-resource pattern, because a static constructor that
    /// eagerly parses several thousand terms would tax every consumer whether or not they resolve
    /// anything.
    /// </summary>
    public sealed class ControlledVocabulary
    {
        private const string PsiMsResourceName = "UsefulProteomicsDatabases.Resources.psi-ms.obo";
        private const string PrideResourceName = "UsefulProteomicsDatabases.Resources.pride_cv.obo";

        private static readonly Lazy<ControlledVocabulary> LazyPsiMs =
            new(() => Load(PsiMsResourceName, "MS"));

        private static readonly Lazy<ControlledVocabulary> LazyPride =
            new(() => Load(PrideResourceName, "PRIDE"));

        private readonly Dictionary<string, CvParam> _byAccession;
        private readonly Dictionary<string, CvParam> _byName;

        private ControlledVocabulary(string label, string version,
            Dictionary<string, CvParam> byAccession, Dictionary<string, CvParam> byName)
        {
            Label = label;
            Version = version;
            _byAccession = byAccession;
            _byName = byName;
        }

        /// <summary>The PSI-MS controlled vocabulary: instruments, dissociation methods, enzymes, analyzers.</summary>
        public static ControlledVocabulary PsiMs => LazyPsiMs.Value;

        /// <summary>
        /// The PRIDE controlled vocabulary. Small but not optional: real annotations use it for
        /// terms PSI-MS does not carry, notably the acquisition method
        /// ("NT=Data-dependent acquisition;AC=PRIDE:0000627").
        /// </summary>
        public static ControlledVocabulary Pride => LazyPride.Value;

        /// <summary>The CV label used in an accession prefix and in a cvParam's cvRef, e.g. "MS".</summary>
        public string Label { get; }

        /// <summary>
        /// The <c>data-version</c> the embedded snapshot was built from, e.g. "4.1.258". Record this
        /// alongside any terms resolved here: without it there is no way to tell whether two
        /// documents disagree because the annotation changed or because the vocabulary did.
        /// Empty if the file declared none.
        /// </summary>
        public string Version { get; }

        /// <summary>
        /// Every term. No ordering is promised — this is a dictionary's value collection, and
        /// <c>Dictionary&lt;K,V&gt;.Values</c> guarantees nothing about order regardless of how it
        /// happens to behave today.
        /// </summary>
        public IReadOnlyCollection<CvParam> Terms => _byAccession.Values;

        /// <summary>
        /// Look a term up by its accession — the stable machine identifier, and the only correct
        /// key. Accession comparison is case-sensitive on the prefix by convention ("MS:1001911"),
        /// but this accepts any casing rather than failing a caller over it.
        /// </summary>
        public bool TryGetByAccession(string accession, out CvParam term)
        {
            term = null;
            return !string.IsNullOrWhiteSpace(accession)
                   && _byAccession.TryGetValue(accession.Trim(), out term);
        }

        /// <summary>
        /// Look a term up by name, including its exact synonyms. Case- and whitespace-insensitive.
        ///
        /// This exists for the one job that genuinely needs it — turning free text a vendor wrote
        /// ("Orbitrap Fusion Lumos" out of a RAW file, "HCD" out of a protocol) into an accession.
        /// It is NOT how terms should be compared: names are display labels and change between
        /// versions, which is why <see cref="MzLibUtil.CvParam"/> documents matching on Accession.
        /// A name that is ambiguous across two terms resolves to the first one in file order, so a
        /// caller that cares should check <see cref="TryGetByAccession"/> against a known id.
        /// </summary>
        public bool TryGetByName(string name, out CvParam term)
        {
            term = null;
            return !string.IsNullOrWhiteSpace(name)
                   && _byName.TryGetValue(Normalize(name), out term);
        }

        private static ControlledVocabulary Load(string resourceName, string label)
        {
            var assembly = typeof(ControlledVocabulary).Assembly;

            using var stream = assembly.GetManifestResourceStream(resourceName)
                ?? throw new MzLibException(
                    $"Could not find embedded resource '{resourceName}'. It should be declared as an " +
                    "EmbeddedResource in UsefulProteomicsDatabases.csproj.");

            // OboParser reads from a path, so the text is materialised and handed to ParseText
            // rather than spilled to a temp file. These files are 1.1 MB and 224 KB.
            using var reader = new StreamReader(stream, Encoding.UTF8);
            string text = reader.ReadToEnd();

            string version = ReadDataVersion(text);

            var byAccession = new Dictionary<string, CvParam>(StringComparer.OrdinalIgnoreCase);
            var byName = new Dictionary<string, CvParam>(StringComparer.Ordinal);

            foreach (var oboTerm in new OboParser().ParseText(text))
            {
                if (string.IsNullOrWhiteSpace(oboTerm.Id) || string.IsNullOrWhiteSpace(oboTerm.Name))
                    continue;

                // Obsolete terms stay resolvable by ACCESSION -- an old document legitimately
                // contains them and must still be readable -- but are kept out of the name index so
                // they are never chosen when authoring something new.
                bool obsolete = oboTerm.ValuePairs?.Any(p =>
                    string.Equals(p.Tag, "is_obsolete", StringComparison.Ordinal)
                    && p.Value.StartsWith("true", StringComparison.OrdinalIgnoreCase)) ?? false;

                // The CV label comes from the ACCESSION, not from which file the term was read out
                // of. An .obo routinely imports terms from other vocabularies: psi-ms.obo carries
                // 128 foreign terms (UO, PEFF, NCIT — e.g. UO:0000266 electronvolt) and
                // pride_cv.obo carries 19, including MS:1000044 itself. Stamping the file's label on
                // all of them would emit cvRef="MS" for a UO term, which is wrong in an mzML
                // cvParam and wrong in an SDRF cell.
                var term = new CvParam(LabelFor(oboTerm.Id, label), oboTerm.Id, oboTerm.Name, "");
                byAccession[oboTerm.Id] = term;

                if (obsolete)
                    continue;

                // First writer wins, so an earlier term keeps an ambiguous name.
                string key = Normalize(oboTerm.Name);
                if (!byName.ContainsKey(key))
                    byName[key] = term;

                foreach (string synonym in ExactSynonyms(oboTerm))
                {
                    string synonymKey = Normalize(synonym);
                    if (!byName.ContainsKey(synonymKey))
                        byName[synonymKey] = term;
                }
            }

            return new ControlledVocabulary(label, version, byAccession, byName);
        }

        /// <summary>
        /// OboTerm has no typed synonym collection: synonyms arrive as untyped tag/value pairs whose
        /// value is the raw OBO line, e.g.
        ///   synonym: "collision-induced dissociation" EXACT []
        /// Only EXACT synonyms are indexed. RELATED and BROAD ones name something adjacent rather
        /// than the same thing, and treating them as equivalent is how a lookup quietly returns a
        /// parent term.
        /// </summary>
        private static IEnumerable<string> ExactSynonyms(OboTerm term)
        {
            if (term.ValuePairs == null)
                yield break;

            foreach (var pair in term.ValuePairs)
            {
                if (!string.Equals(pair.Tag, "synonym", StringComparison.Ordinal))
                    continue;

                int open = pair.Value.IndexOf('"');
                if (open < 0) continue;
                int close = pair.Value.IndexOf('"', open + 1);
                if (close <= open) continue;

                if (pair.Value.IndexOf("EXACT", close, StringComparison.Ordinal) < 0)
                    continue;

                string synonym = pair.Value.Substring(open + 1, close - open - 1);
                if (!string.IsNullOrWhiteSpace(synonym))
                    yield return synonym;
            }
        }

        /// <summary>
        /// The vocabulary an accession belongs to, taken from its prefix ("UO:0000266" -> "UO").
        /// Falls back to the containing file's label for a malformed id.
        /// </summary>
        private static string LabelFor(string accession, string fileLabel)
        {
            int colon = accession.IndexOf(':');
            return colon > 0 ? accession.Substring(0, colon) : fileLabel;
        }

        private static string ReadDataVersion(string text)
        {
            // Scan only the header. Splitting the whole file materialised ~29,500 strings to read
            // what is line 2, and the break saved nothing because Split had already done the work.
            const string tag = "data-version:";
            int start = 0;
            while (start < text.Length)
            {
                int end = text.IndexOf('\n', start);
                if (end < 0) end = text.Length;

                if (end - start >= tag.Length
                    && string.CompareOrdinal(text, start, tag, 0, tag.Length) == 0)
                    return text.Substring(start + tag.Length, end - start - tag.Length).Trim();

                // The first stanza header ends the file header. Both files open with [Typedef]
                // rather than [Term], so stop at either.
                if (text[start] == '[')
                    break;

                start = end + 1;
            }
            return "";
        }

        /// <summary>
        /// Case- and whitespace-insensitive key for name lookup. Conservative on purpose: it
        /// collapses runs of whitespace and lowercases, and does nothing else, so it cannot merge
        /// two genuinely different terms.
        /// </summary>
        private static string Normalize(string value) =>
            string.Join(' ', value.Split((char[])null, StringSplitOptions.RemoveEmptyEntries))
                .ToLowerInvariant();
    }
}
