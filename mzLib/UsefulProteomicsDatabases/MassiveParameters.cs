using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml;
using MzLibUtil;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// What mzLib reads from a MassIVE dataset's <c>ccms_parameters/params.xml</c>, the record of its submission form:
    /// the map from MassIVE's stored file names to the paths the submitter uploaded, and the species and instruments
    /// the submitter chose.
    /// </summary>
    /// <remarks>
    /// <b>An allow-list, never more.</b> The file also holds the submitter's email address, the PI's name, email and
    /// institution, and, on older datasets, a reviewer password (<c>dataset.password</c>, MSV000080763). None of those
    /// is read into any property, so none can leave this parser.
    /// <para>
    /// It holds no protease, tolerance or label (measured 2026-10-18 on MSV000096226), so a search's settings must
    /// come from elsewhere. An update has a <c>params.xml</c> of its own (<c>updates/&lt;name&gt;/ccms_parameters/</c>),
    /// with <see cref="Workflow"/> "UPDATE-MASSIVE-DATASET"; it maps only the update's own files.
    /// </para>
    /// <para>
    /// The null rule is <see cref="MassiveDataset"/>'s: a property is null when the file has no such parameter, and
    /// empty when the parameter is there with nothing in it.
    /// </para>
    /// </remarks>
    public sealed class MassiveParameters
    {
        private MassiveParameters()
        {
        }

        /// <summary>The submission's workflow: "MASSIVE-COMPLETE" for a dataset, "UPDATE-MASSIVE-DATASET" for an update; null if absent.</summary>
        public string Workflow { get; private init; }

        /// <summary>
        /// Every <c>upload_file_mapping</c>, in file order: MassIVE's stored name ("RAW-00000.raw") and the path the
        /// submitter uploaded it from. Null when the file has none.
        /// </summary>
        public IReadOnlyList<MassiveUploadMapping> UploadFileMappings { get; private init; }

        /// <summary>The species the submitter chose (<c>dataset.species</c>), e.g. "NCBITaxon:9606"; null when absent.</summary>
        public IReadOnlyList<string> Species { get; private init; }

        /// <summary>The instruments the submitter chose (<c>dataset.instrument</c>), e.g. "MS:1003028"; null when absent.</summary>
        public IReadOnlyList<string> Instruments { get; private init; }

        /// <summary>Reads a <c>params.xml</c> file.</summary>
        /// <exception cref="ArgumentException">The path is blank.</exception>
        /// <exception cref="MzLibException">The file is not a ProteoSAFe parameters file, or an upload mapping has no '|'.</exception>
        public static MassiveParameters Read(string path)
        {
            if (string.IsNullOrWhiteSpace(path))
                throw new ArgumentException("A params.xml path is required.", nameof(path));
            using FileStream stream = File.OpenRead(path);
            return Read(stream);
        }

        /// <summary>Reads a <c>params.xml</c> from a stream. DTDs are refused and no external entity is resolved.</summary>
        /// <exception cref="ArgumentNullException">The stream is null.</exception>
        /// <exception cref="MzLibException">The stream is not a ProteoSAFe parameters file, or an upload mapping has no '|'.</exception>
        public static MassiveParameters Read(Stream stream)
        {
            if (stream == null)
                throw new ArgumentNullException(nameof(stream));

            var settings = new XmlReaderSettings { DtdProcessing = DtdProcessing.Prohibit, XmlResolver = null, IgnoreComments = true };
            var values = new List<(string Name, string Value)>();
            try
            {
                using XmlReader reader = XmlReader.Create(stream, settings);
                if (reader.MoveToContent() != XmlNodeType.Element || reader.Name != "parameters")
                    throw new MzLibException($"This is not a ProteoSAFe params.xml: its root element is <{reader.Name}>, not <parameters>.");

                reader.Read();
                while (!reader.EOF)
                {
                    if (reader.NodeType == XmlNodeType.Element && reader.Name == "parameter"
                        && reader.GetAttribute("name") is "upload_file_mapping" or "workflow" or "dataset.species" or "dataset.instrument")
                    {
                        string name = reader.GetAttribute("name");
                        // Reading the content moves the reader past the element already: no Read() after it, or the
                        // next <parameter> (the next mapping, say) would be skipped.
                        values.Add((name, reader.ReadElementContentAsString()));
                        continue;
                    }
                    reader.Read();
                }
            }
            catch (XmlException e)
            {
                throw new MzLibException("This is not a readable ProteoSAFe params.xml.", e);
            }

            List<string> Of(string name) => values.Where(v => v.Name == name).Select(v => v.Value).ToList();
            List<string> mappings = Of("upload_file_mapping");
            List<string> species = Of("dataset.species");
            List<string> instruments = Of("dataset.instrument");

            return new MassiveParameters
            {
                Workflow = Of("workflow").FirstOrDefault(),
                UploadFileMappings = mappings.Count == 0 ? null : mappings.Select(ParseMapping).ToList(),
                Species = species.Count == 0 ? null : SplitList(species),
                Instruments = instruments.Count == 0 ? null : SplitList(instruments),
            };
        }

        /// <summary>"RAW-00000.raw|user/folder/run.raw": the stored name, then the uploaded path.</summary>
        private static MassiveUploadMapping ParseMapping(string value)
        {
            int bar = value.IndexOf('|');
            if (bar <= 0)
                throw new MzLibException($"An upload_file_mapping in params.xml has no stored name before a '|': '{value}'.");
            return new MassiveUploadMapping(value.Substring(0, bar), value.Substring(bar + 1));
        }

        /// <summary>The values of ';'-separated parameters ("NCBITaxon:9606;"), trimmed, empty ones dropped.</summary>
        private static List<string> SplitList(IEnumerable<string> values) =>
            values.SelectMany(v => v.Split(';')).Select(v => v.Trim()).Where(v => v.Length > 0).ToList();
    }

    /// <summary>One <c>upload_file_mapping</c> of a MassIVE <c>params.xml</c>.</summary>
    /// <param name="StoredName">MassIVE's name for the file, e.g. "RAW-00000.raw" or "metadata-00000.tsv".</param>
    /// <param name="UploadPath">The path the submitter uploaded it from, e.g. "ccms/dataset/pride_imports/PXD004319/RAW/x.raw".</param>
    public sealed record MassiveUploadMapping(string StoredName, string UploadPath);
}
