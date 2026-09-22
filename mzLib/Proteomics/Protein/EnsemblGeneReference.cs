using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace Proteomics
{
    /// <summary>
    /// One Ensembl transcript a protein entry links to, and the gene that transcript belongs to.
    ///
    /// This is a derived view over Protein.DatabaseReferences, not a stored field -- UniProt's
    /// &lt;dbReference type="Ensembl"&gt; entries already arrive there because ProteinXmlEntry keeps
    /// every dbReference generically. See Protein.EnsemblGeneReferences and Protein.EnsemblGeneIds.
    ///
    /// UniProt writes one reference per transcript: the reference's Id is the transcript
    /// ("ENST00000229239.10"), and the gene is a property of it
    /// (&lt;property type="gene ID" value="ENSG00000111640.15"/&gt;). That is why this is not a
    /// FirstOrDefault(...)?.Id projection like Protein.NcbiTaxonomyId: that shape would return a
    /// transcript, and would drop every gene after the first.
    /// </summary>
    public class EnsemblGeneReference
    {
        /// <summary>The property type UniProt uses for the gene id on an Ensembl dbReference.</summary>
        public const string GeneIdPropertyType = "gene ID";

        /// <summary>The property type UniProt uses for the protein id on an Ensembl dbReference.</summary>
        public const string ProteinIdPropertyType = "protein sequence ID";

        public EnsemblGeneReference(string transcriptId, string proteinId, string versionedGeneId)
        {
            TranscriptId = transcriptId ?? "";
            ProteinId = proteinId ?? "";
            VersionedGeneId = versionedGeneId ?? "";
            SplitVersion(VersionedGeneId, out string geneId, out int? version);
            GeneId = geneId;
            GeneVersion = version;
        }

        /// <summary>The Ensembl transcript id as UniProt wrote it, versioned, e.g. "ENST00000229239.10".</summary>
        public string TranscriptId { get; }

        /// <summary>The Ensembl protein id as UniProt wrote it, versioned. Empty when absent.</summary>
        public string ProteinId { get; }

        /// <summary>
        /// The gene id exactly as UniProt wrote it, e.g. "ENSG00000111640.15". The only form that pins
        /// the link to a specific Ensembl release, so it is kept for audit.
        /// </summary>
        public string VersionedGeneId { get; }

        /// <summary>
        /// The stable gene id, e.g. "ENSG00000111640" -- the one to join on, because versions change for
        /// reasons that are not biology. Equal to VersionedGeneId when there is no numeric version.
        /// </summary>
        public string GeneId { get; }

        /// <summary>The numeric gene version, or null when the id carried none. Absent is not version 0.</summary>
        public int? GeneVersion { get; }

        /// <summary>
        /// Projects a flat list of database references onto Ensembl gene links: keeps only
        /// Type == "Ensembl" references that carry a gene id, and keeps each transcript once.
        ///
        /// Filtering by type is mandatory, not defensive -- the same list holds GO, PubMed, RefSeq and
        /// organism references, and the XML parser applies no depth guard. A reference without a gene
        /// id is not a gene link and is left out rather than reported with an empty gene.
        /// </summary>
        public static IReadOnlyList<EnsemblGeneReference> FromDatabaseReferences(IEnumerable<DatabaseReference> references)
        {
            if (references == null)
            {
                return Array.Empty<EnsemblGeneReference>();
            }

            var seenTranscripts = new HashSet<string>(StringComparer.Ordinal);
            var result = new List<EnsemblGeneReference>();
            foreach (var reference in references)
            {
                if (reference == null || reference.Type != Protein.EnsemblDatabaseReferenceType)
                {
                    continue;
                }

                // Matched by type, never by position: ProteinDbWriter re-sorts properties on write.
                string geneId = PropertyValue(reference, GeneIdPropertyType);
                if (string.IsNullOrEmpty(geneId) || !seenTranscripts.Add(reference.Id ?? ""))
                {
                    continue;
                }

                result.Add(new EnsemblGeneReference(reference.Id, PropertyValue(reference, ProteinIdPropertyType), geneId));
            }

            return result;
        }

        private static string PropertyValue(DatabaseReference reference, string type) =>
            (reference.Properties ?? Enumerable.Empty<Tuple<string, string>>())
                .FirstOrDefault(p => p != null && p.Item1 == type)?.Item2;

        /// <summary>
        /// Splits "ENSG00000111640.15" into "ENSG00000111640" and 15. Only an all-digit suffix after the
        /// last '.' is a version; anything else is left whole with no version rather than guessed at.
        /// </summary>
        private static void SplitVersion(string versionedId, out string stableId, out int? version)
        {
            stableId = versionedId;
            version = null;

            int dot = versionedId.LastIndexOf('.');
            if (dot <= 0 || dot == versionedId.Length - 1)
            {
                return;
            }

            string suffix = versionedId.Substring(dot + 1);
            if (suffix.All(char.IsAsciiDigit)
                && int.TryParse(suffix, NumberStyles.None, CultureInfo.InvariantCulture, out int parsed))
            {
                stableId = versionedId.Substring(0, dot);
                version = parsed;
            }
        }
    }
}
