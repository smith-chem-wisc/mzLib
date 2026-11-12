using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using MzLibUtil;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// A MassIVE dataset's record, as MassIVE's PROXI service gives it (<see cref="MassiveClient.GetDatasetAsync"/>).
    /// </summary>
    /// <remarks>
    /// <b>One rule for every field.</b> Null means the record does not state it: the field is missing, its value is
    /// JSON null, or its value is the string "null" (MSV000078494 says "null" for its species and its instrument).
    /// Empty (an empty string, or a list with no items) means the record states it as empty. A dataset MassIVE does
    /// not know has no record at all: <see cref="MassiveClient.TryGetDatasetAsync"/> answers <c>Found = false</c>.
    /// So "not stated" and "stated as nothing" stay apart, as aging (029) and sdrf (026) asked.
    /// <para>
    /// Contacts are deliberately not exposed: no consumer asked for them, and they carry email addresses.
    /// </para>
    /// </remarks>
    public sealed class MassiveDataset
    {
        /// <summary>The PSI-MS accession of a MassIVE dataset identifier term.</summary>
        internal const string MassiveIdentifierAccession = "MS:1002487";

        /// <summary>The PSI-MS accession of a ProteomeXchange accession number term.</summary>
        internal const string ProteomeXchangeAccession = "MS:1001919";

        /// <summary>The MassIVE accession, e.g. "MSV000096226". Never null.</summary>
        public string MassiveAccession { get; init; }

        /// <summary>The ProteomeXchange accession, e.g. "PXD057243", or null when the dataset has none.</summary>
        public string PxdAccession { get; init; }

        /// <summary>
        /// The title, with a leading "Title: " cut (see <see cref="TitleHadPrefix"/>), or null when not stated.
        /// </summary>
        public string Title { get; init; }

        /// <summary>
        /// Whether the record's title began "Title: ", which <see cref="Title"/> has cut. It is a submitter's habit,
        /// not the server's (1 title of 121 measured, MSV000097834), and it must never reach a drafted SDRF (sdrf MV2-T4).
        /// </summary>
        public bool TitleHadPrefix { get; init; }

        /// <summary>The summary, verbatim, or null when not stated.</summary>
        public string Summary { get; init; }

        /// <summary>The submitter keywords (MS:1001925 values), in the record's order, or null when not stated.</summary>
        public IReadOnlyList<string> Keywords { get; init; }

        /// <summary>
        /// The species, one entry per term group, or null when not stated. A group whose every value is unstated is
        /// left out, so MSV000078494 (whose only species term says "null") gives an empty list, not null.
        /// </summary>
        public IReadOnlyList<MassiveSpecies> Species { get; init; }

        /// <summary>
        /// The instruments, or null when not stated. The generic "instrument model" term (MS:1000031) is kept when it
        /// names one in its value ("Bruker timsTOF Pro": 7 of the 8 such terms in 121 records do) and left out when
        /// its value is unstated.
        /// </summary>
        public IReadOnlyList<CvParam> Instruments { get; init; }

        /// <summary>The modifications (UNIMOD terms), or null when not stated.</summary>
        public IReadOnlyList<CvParam> Modifications { get; init; }

        /// <summary>
        /// The publications, one entry per term group, or null when not stated. A record that says it has no
        /// associated manuscript (MS:1002853) gives an empty list.
        /// </summary>
        public IReadOnlyList<MassivePublication> Publications { get; init; }

        /// <summary>The dataset's MassIVE web page (MS:1002488), or null when not stated.</summary>
        public string DatasetUri { get; init; }

        /// <summary>
        /// The dataset's folder on MassIVE's FTP server, from the record's FTP location (MS:1002852), e.g.
        /// "v07/MSV000096226", or null when not stated. Only the path is kept: the record names the host
        /// <c>massive.ucsd.edu</c>, which does not answer FTP, and the files are on <see cref="MassiveClient.FtpHost"/>.
        /// </summary>
        public string FtpLocation { get; init; }

        // ---- PROXI parsing ------------------------------------------------------------------------------------

        private const string KeywordAccession = "MS:1001925";
        private const string TaxIdAccession = "MS:1001467";
        private const string CommonNameAccession = "MS:1001468";
        private const string ScientificNameAccession = "MS:1001469";
        private const string InstrumentModelAccession = "MS:1000031";
        private const string PubMedAccession = "MS:1000879";
        private const string DoiAccession = "MS:1001922";
        private const string ReferenceAccession = "MS:1002866";
        private const string DatasetUriAccession = "MS:1002488";
        private const string FtpLocationAccession = "MS:1002852";
        private const string TitlePrefix = "Title: ";

        private static readonly Encoding StrictUtf8 = new UTF8Encoding(encoderShouldEmitUTF8Identifier: false, throwOnInvalidBytes: true);

        private static readonly Lazy<Encoding> Windows1252 = new(() =>
        {
            // Inbox on .NET (System.Text.Encoding.CodePages ships with the shared framework); registering twice is harmless.
            Encoding.RegisterProvider(CodePagesEncodingProvider.Instance);
            return Encoding.GetEncoding(1252);
        });

        /// <summary>
        /// Decodes a PROXI body: strict UTF-8, else windows-1252. The server's charset label is ignored, because every
        /// record says ISO-8859-1 while the bytes are mixed: MSV000088323 is Latin-1 ("Rapp\xE9"), MSV000079927 is
        /// UTF-8 ("CAUTI\xC2\xB4s"), and decoding the second by its label turns its one acute accent into two
        /// characters (measured 2026-11-03).
        /// windows-1252, not Latin-1: the two differ at 0x80-0x9F, where 1252 has curly quotes and dashes.
        /// </summary>
        internal static string DecodeProxiBody(byte[] body)
        {
            string text;
            try
            {
                text = StrictUtf8.GetString(body);
            }
            catch (DecoderFallbackException)
            {
                text = Windows1252.Value.GetString(body);
            }
            return text.Length > 0 && text[0] == '﻿' ? text.Substring(1) : text;
        }

        /// <summary>
        /// Reads a MassIVE PROXI dataset body. A body that is not a JSON object, or carries no MassIVE identifier, is
        /// not a record, and throws <see cref="MzLibException"/>.
        /// </summary>
        /// <param name="body">The response bytes.</param>
        /// <param name="requested">The accession asked for, for the message.</param>
        internal static MassiveDataset FromProxi(byte[] body, string requested)
        {
            JObject record;
            try
            {
                using var reader = new JsonTextReader(new StringReader(DecodeProxiBody(body)))
                {
                    // Values are kept as the strings MassIVE sent: a date-shaped one must not come back reformatted.
                    DateParseHandling = DateParseHandling.None,
                    FloatParseHandling = FloatParseHandling.Decimal,
                };
                record = JToken.ReadFrom(reader) as JObject;
            }
            catch (JsonReaderException e)
            {
                throw new MzLibException($"MassIVE's PROXI record for '{requested}' is not valid JSON.", e);
            }

            if (record == null)
                throw new MzLibException($"MassIVE's PROXI answer for '{requested}' is not a dataset record.");

            List<JObject> identifiers = Terms(record["accession"]).ToList();
            string massiveAccession = ValueOf(identifiers, MassiveIdentifierAccession);
            if (massiveAccession == null)
                throw new MzLibException($"MassIVE's PROXI record for '{requested}' carries no MassIVE identifier ({MassiveIdentifierAccession}).");

            string title = Stated(record["title"]);
            bool titleHadPrefix = title != null && title.StartsWith(TitlePrefix, StringComparison.Ordinal);
            if (titleHadPrefix)
                title = title.Substring(TitlePrefix.Length);

            List<JObject> links = Terms(record["datasetLink"]).ToList();

            return new MassiveDataset
            {
                MassiveAccession = massiveAccession,
                PxdAccession = ValueOf(identifiers, ProteomeXchangeAccession),
                Title = title,
                TitleHadPrefix = titleHadPrefix,
                Summary = Stated(record["summary"]),
                Keywords = ListOf(record["keywords"], field => Terms(field)
                    .Where(t => Accession(t) == KeywordAccession)
                    .Select(t => Stated(t["value"]))
                    .Where(v => v != null)),
                Species = ListOf(record["species"], field => Groups(field)
                    .Select(g => new MassiveSpecies(ValueOf(g, TaxIdAccession), ValueOf(g, ScientificNameAccession), ValueOf(g, CommonNameAccession)))
                    .Where(s => s.TaxId != null || s.ScientificName != null || s.CommonName != null)),
                Instruments = ListOf(record["instruments"], field => Terms(field)
                    .Where(t => Accession(t) != InstrumentModelAccession || Stated(t["value"]) != null)
                    .Select(ToCvParam)),
                Modifications = ListOf(record["modifications"], field => Terms(field).Select(ToCvParam)),
                Publications = ListOf(record["publications"], field => Groups(field)
                    .Select(g => new MassivePublication(ValueOf(g, PubMedAccession), ValueOf(g, DoiAccession), ValueOf(g, ReferenceAccession)))
                    .Where(p => p.PubMedId != null || p.Doi != null || p.Reference != null)),
                DatasetUri = ValueOf(links, DatasetUriAccession),
                FtpLocation = FtpPathOf(ValueOf(links, FtpLocationAccession)),
            };
        }

        /// <summary>
        /// A string value under the null rule: null for a missing token, JSON null, or the string "null"; otherwise
        /// the value as sent, an empty string included.
        /// </summary>
        private static string Stated(JToken token)
        {
            if (token == null || token.Type is JTokenType.Null or JTokenType.Undefined)
                return null;
            string value = token.Type == JTokenType.String ? (string)token : token.ToString(Formatting.None);
            return string.Equals(value.Trim(), "null", StringComparison.OrdinalIgnoreCase) ? null : value;
        }

        /// <summary>A list field under the null rule: null when not stated, else the items <paramref name="read"/> gives.</summary>
        private static IReadOnlyList<T> ListOf<T>(JToken field, Func<JToken, IEnumerable<T>> read) =>
            field == null || field.Type is JTokenType.Null or JTokenType.Undefined ? null : read(field).ToList();

        /// <summary>The terms of a field that holds one term or a list of terms; groups inside it are flattened.</summary>
        private static IEnumerable<JObject> Terms(JToken field) => field switch
        {
            JObject term => new[] { term },
            JArray items => items.SelectMany(Terms),
            _ => Enumerable.Empty<JObject>(),
        };

        /// <summary>
        /// The term groups of a field. PROXI gives a list of groups (each a list of terms), but a record can also give
        /// bare terms (71 of 121 records give their publications as one bare term, "no associated manuscript"); each
        /// bare term is then a group of its own.
        /// </summary>
        private static IEnumerable<List<JObject>> Groups(JToken field) => field switch
        {
            JObject term => new[] { new List<JObject> { term } },
            JArray items => items.Select(item => Terms(item).ToList()).Where(g => g.Count > 0),
            _ => Enumerable.Empty<List<JObject>>(),
        };

        private static string Accession(JObject term) => (string)term["accession"];

        /// <summary>The stated value of the first term in <paramref name="terms"/> with <paramref name="accession"/>, else null.</summary>
        private static string ValueOf(IEnumerable<JObject> terms, string accession) =>
            terms.Where(t => Accession(t) == accession).Select(t => Stated(t["value"])).FirstOrDefault(v => v != null);

        private static CvParam ToCvParam(JObject term) => new(
            (string)term["cvLabel"] ?? string.Empty,
            Accession(term) ?? string.Empty,
            (string)term["name"] ?? string.Empty,
            Stated(term["value"]) ?? string.Empty);

        /// <summary>
        /// The path part of the record's FTP location, without its leading '/', or null. Taken from the text, not
        /// through <see cref="Uri"/>, which would resolve a ".." away before anything could refuse it.
        /// </summary>
        private static string FtpPathOf(string ftpLocation)
        {
            if (ftpLocation == null)
                return null;
            int scheme = ftpLocation.IndexOf("://", StringComparison.Ordinal);
            string path = ftpLocation;
            if (scheme >= 0)
            {
                int slash = ftpLocation.IndexOf('/', scheme + 3);
                path = slash >= 0 ? ftpLocation.Substring(slash) : string.Empty;
            }
            return Uri.UnescapeDataString(path).Trim('/');
        }
    }

    /// <summary>One species of a MassIVE dataset: one PROXI term group. Each value is null when the group does not state it.</summary>
    /// <param name="TaxId">The NCBI taxonomy ID (MS:1001467), e.g. "10090".</param>
    /// <param name="ScientificName">The scientific name (MS:1001469), e.g. "Mus musculus".</param>
    /// <param name="CommonName">The common name (MS:1001468).</param>
    public sealed record MassiveSpecies(string TaxId, string ScientificName, string CommonName);

    /// <summary>One publication of a MassIVE dataset: one PROXI term group. Each value is null when the group does not state it.</summary>
    /// <param name="PubMedId">The PubMed identifier (MS:1000879).</param>
    /// <param name="Doi">The DOI (MS:1001922).</param>
    /// <param name="Reference">The reference text (MS:1002866).</param>
    public sealed record MassivePublication(string PubMedId, string Doi, string Reference);
}
