using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Xml;
using Easy.Common.Extensions;
using Omics;
using Omics.BioPolymer;
using Omics.Modifications;
using Transcriptomics;
using System.Data;

namespace UsefulProteomicsDatabases
{

    /// <summary>
    /// Provides methods for writing protein and nucleic acid databases to XML and FASTA formats.
    /// Did not rename to DbWriter to ensure compatibility with the original UsefulProteomicsDatabases namespace.
    /// </summary>
    public class ProteinDbWriter
    {
        /// <summary>
        /// Writes an XML database for a list of RNA sequences, including additional modifications.
        /// </summary>
        /// <param name="additionalModsToAddToProteins">A dictionary of additional modifications to add to proteins.</param>
        /// <param name="bioPolymerList">A list of RNA sequences to be written to the database.</param>
        /// <param name="outputFileName">The name of the output XML file.</param>
        /// <returns>A dictionary of new modification residue entries.</returns>
        public static Dictionary<string, int> WriteXmlDatabase(
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins,
            List<IBioPolymer> bioPolymerList,
            string outputFileName)
        {
            return bioPolymerList.Any(p => p is Protein)
                ? WriteXmlDatabase(additionalModsToAddToProteins, bioPolymerList.Cast<Protein>().ToList(), outputFileName)
                : WriteXmlDatabase(additionalModsToAddToProteins, bioPolymerList.Cast<RNA>().ToList(), outputFileName);
        }

        /// <summary>
        /// Writes an XML database for a list of nucleic acid sequences, including additional modifications.
        /// </summary>
        /// <param name="additionalModsToAddToNucleicAcids">A dictionary of additional modifications to add to proteins.</param>
        /// <param name="nucleicAcidList">A list of nucleic acid sequences to be written to the database.</param>
        /// <param name="outputFileName">The name of the output XML file.</param>
        /// <param name="updateTimeStamp">If true, updates the modified attribute to today's date when attributes are written (currently RNA omits attributes as per original).</param>
        /// <returns>A dictionary of new modification residue entries.</returns>
        /// <remarks>
        /// Several chunks of code are commented out. These are blocks that are intended to be implemented in the future, but
        /// are not necessary for the bare bones implementation of Transcriptomics
        /// </remarks>
        public static Dictionary<string, int> WriteXmlDatabase(
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToNucleicAcids,
            List<RNA> nucleicAcidList,
            string outputFileName,
            bool updateTimeStamp = false)
        {
            additionalModsToAddToNucleicAcids ??= new Dictionary<string, HashSet<Tuple<int, Modification>>>();

            // Write non-variant RNA (when variants aren't applied, this just returns the RNA itself)
            var nonVariantRna = nucleicAcidList.Select(p => p.ConsensusVariant).OfType<RNA>().Distinct().ToList();

            Dictionary<string, int> newModResEntries = new();

            using (XmlWriter writer = XmlWriter.Create(outputFileName, CreateIndentedWriterSettings()))
            {
                WriteStartDocument(writer);

                // Modifications catalog
                var allRelevantMods = CollectAllRelevantModsForRna(nonVariantRna, additionalModsToAddToNucleicAcids);
                WriteModificationCatalog(writer, allRelevantMods);

                // Entries
                foreach (var rna in nonVariantRna)
                {
                    WriteRnaEntry(writer, rna, additionalModsToAddToNucleicAcids, newModResEntries, updateTimeStamp);
                }

                WriteEndDocument(writer);
            }

            return newModResEntries;
        }

        /// <summary>
        /// Writes a protein database in mzLibProteinDb format, with additional modifications from the AdditionalModsToAddToProteins list.
        /// </summary>
        /// <param name="additionalModsToAddToProteins"></param>
        /// <param name="proteinList"></param>
        /// <param name="outputFileName"></param>
        /// <param name="updateTimeStamp"></param>
        /// <param name="includeAppliedVariantEntries">
        /// If true, applied (realized) variant proteoforms (with a different accession produced by VariantApplication) are written
        /// as separate &lt;entry&gt; elements in addition to their consensus (canonical) parents.
        /// </param>
        /// <param name="includeAppliedVariantFeatures">
        /// If true and an applied variant entry is written, its AppliedSequenceVariations are emitted as
        /// &lt;feature type="sequence variant"&gt; elements so differences remain explicit (even though its BaseSequence already contains them).
        /// </param>
        /// <returns>The new "modified residue" entries that are added due to being in the Mods dictionary</returns>
        public static Dictionary<string, int> WriteXmlDatabase(
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins,
            List<Protein> proteinList,
            string outputFileName,
            bool updateTimeStamp = false,
            bool includeAppliedVariantEntries = false,
            bool includeAppliedVariantFeatures = true)
        {
            additionalModsToAddToProteins ??= new Dictionary<string, HashSet<Tuple<int, Modification>>>();

            var proteinsToWrite = BuildProteinsToWrite(proteinList, includeAppliedVariantEntries);

            Dictionary<string, int> newModResEntries = new();

            using (XmlWriter writer = XmlWriter.Create(outputFileName, CreateIndentedWriterSettings()))
            {
                WriteStartDocument(writer);

                // Modifications catalog
                var allRelevantMods = CollectAllRelevantModsForProteins(proteinsToWrite, includeAppliedVariantEntries, additionalModsToAddToProteins);
                WriteModificationCatalog(writer, allRelevantMods);

                // Entries
                foreach (var protein in proteinsToWrite.OrderBy(p => p.Accession, StringComparer.Ordinal))
                {
                    bool isAppliedVariantEntry = DetermineIsAppliedVariantEntry(protein, includeAppliedVariantEntries);
                    WriteProteinEntry(writer, protein, isAppliedVariantEntry, updateTimeStamp, includeAppliedVariantFeatures, additionalModsToAddToProteins, newModResEntries);
                }

                WriteEndDocument(writer);
            }

            return newModResEntries;
        }

        /// <summary>
        /// Writes a FASTA file for a list of proteins.
        /// </summary>
        public static void WriteFastaDatabase(List<Protein> proteinList, string outputFileName, string delimeter)
        {
            using (StreamWriter writer = new StreamWriter(outputFileName))
            {
                foreach (Protein protein in proteinList)
                {
                    string header = delimeter == " " ? protein.GetEnsemblFastaHeader() : protein.GetUniProtFastaHeader();
                    writer.WriteLine(">" + header);
                    writer.WriteLine(protein.BaseSequence);
                }
            }
        }

        /// <summary>
        /// Collects all relevant modifications for RNA: base mods, sequence-variant mods, and additional mods scoped by accession keys.
        /// </summary>
        private static IEnumerable<Modification> CollectAllRelevantModsForRna(
            List<RNA> nonVariantRna,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToNucleicAcids)
        {
            HashSet<Modification> allRelevant = new();

            foreach (var p in nonVariantRna)
            {
                // Variant-specific mods
                if (p.SequenceVariations != null)
                {
                    foreach (var sv in p.SequenceVariations)
                    {
                        if (sv?.OneBasedModifications == null) continue;
                        foreach (var kv in sv.OneBasedModifications)
                        {
                            if (kv.Value == null) continue;
                            foreach (var m in kv.Value)
                            {
                                if (m != null) allRelevant.Add(m);
                            }
                        }
                    }
                }

                // Base possible localized mods
                if (p.OneBasedPossibleLocalizedModifications != null)
                {
                    foreach (var kv in p.OneBasedPossibleLocalizedModifications)
                    {
                        if (kv.Value == null) continue;
                        foreach (var m in kv.Value)
                        {
                            if (m != null) allRelevant.Add(m);
                        }
                    }
                }
            }

            // Additional externally supplied mods (keys that match base accession or variant-accession)
            var allowedAccessions = new HashSet<string>(
                nonVariantRna.SelectMany(p =>
                    (p.SequenceVariations ?? new List<SequenceVariation>())
                        .Select(sv => VariantApplication.GetAccession(p, new[] { sv }))
                        .Concat(new[] { p.Accession })),
                StringComparer.Ordinal);

            foreach (var kv in (additionalModsToAddToNucleicAcids ?? new()).Where(kv => allowedAccessions.Contains(kv.Key)))
            {
                foreach (var t in kv.Value)
                {
                    if (t?.Item2 != null) allRelevant.Add(t.Item2);
                }
            }

            return allRelevant.OrderBy(m => m.IdWithMotif);
        }

        /// <summary>
        /// Collects all relevant modifications for proteins: base mods, sequence-variant mods, applied-variant mods (optional), and additional mods by accession.
        /// </summary>
        private static IEnumerable<Modification> CollectAllRelevantModsForProteins(
            List<Protein> proteinsToWrite,
            bool includeAppliedVariantEntries,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins)
        {
            HashSet<Modification> allRelevantModifications = new();

            foreach (var prot in proteinsToWrite)
            {
                if (prot == null) continue;

                // Base possible localized mods
                if (prot.OneBasedPossibleLocalizedModifications != null)
                {
                    foreach (var kv in prot.OneBasedPossibleLocalizedModifications)
                    {
                        if (kv.Value == null) continue;
                        foreach (var m in kv.Value)
                        {
                            if (m != null) allRelevantModifications.Add(m);
                        }
                    }
                }

                // Candidate sequence variants
                if (prot.SequenceVariations != null)
                {
                    foreach (var sv in prot.SequenceVariations)
                    {
                        if (sv?.OneBasedModifications == null) continue;
                        foreach (var kv in sv.OneBasedModifications)
                        {
                            if (kv.Value == null) continue;
                            foreach (var m in kv.Value)
                            {
                                if (m != null) allRelevantModifications.Add(m);
                            }
                        }
                    }
                }

                // Applied sequence variants (when writing applied variant entries)
                if (includeAppliedVariantEntries && prot.AppliedSequenceVariations != null)
                {
                    foreach (var sv in prot.AppliedSequenceVariations)
                    {
                        if (sv?.OneBasedModifications == null) continue;
                        foreach (var kv in sv.OneBasedModifications)
                        {
                            if (kv.Value == null) continue;
                            foreach (var m in kv.Value)
                            {
                                if (m != null) allRelevantModifications.Add(m);
                            }
                        }
                    }
                }
            }

            // Additional externally supplied mods (filter by accession we actually write)
            var accessionsToWrite = new HashSet<string>(proteinsToWrite.Select(p => p.Accession), StringComparer.Ordinal);
            foreach (var kv in additionalModsToAddToProteins.Where(kv => accessionsToWrite.Contains(kv.Key)))
            {
                foreach (var tup in kv.Value)
                {
                    if (tup?.Item2 != null) allRelevantModifications.Add(tup.Item2);
                }
            }

            return allRelevantModifications.OrderBy(m => m.IdWithMotif);
        }

        /// <summary>
        /// Writes the global catalog of modifications required for all entries in the file.
        /// </summary>
        private static void WriteModificationCatalog(XmlWriter writer, IEnumerable<Modification> modifications)
        {
            foreach (Modification mod in modifications)
            {
                writer.WriteStartElement("modification");
                writer.WriteString(mod.ToString() + Environment.NewLine + "//");
                writer.WriteEndElement();
            }
        }

        /// <summary>
        /// Builds the list of proteins to write: canonical consensus entries plus optional applied variant proteoforms.
        /// </summary>
        private static List<Protein> BuildProteinsToWrite(IEnumerable<Protein> proteinList, bool includeAppliedVariantEntries)
        {
            var consensusProteins = proteinList
                .Select(p => p?.ConsensusVariant)
                .OfType<Protein>()
                .Distinct()
                .ToList();

            List<Protein> proteinsToWrite = new(consensusProteins);

            if (!includeAppliedVariantEntries)
            {
                return proteinsToWrite;
            }

            foreach (var p in proteinList)
            {
                if (p == null) continue;
                var consensus = p.ConsensusVariant as Protein;

                bool isAppliedVariant = p.AppliedSequenceVariations != null
                                        && p.AppliedSequenceVariations.Count > 0
                                        && (consensus == null || !ReferenceEquals(p, consensus));

                if (isAppliedVariant && !proteinsToWrite.Any(x => string.Equals(x.Accession, p.Accession, StringComparison.Ordinal)))
                {
                    proteinsToWrite.Add(p);
                }
            }

            return proteinsToWrite;
        }

        /// <summary>
        /// Writes a complete RNA entry (accession, names, gene/organism, features, sequence).
        /// </summary>
        private static void WriteRnaEntry(
            XmlWriter writer,
            RNA rna,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalMods,
            Dictionary<string, int> newModResEntries,
            bool updateTimeStamp)
        {
            writer.WriteStartElement("entry", "undefined"); // placeholder to match original behavior

            // Accession
            WriteAccession(writer, rna.Accession);

            // Optional presentation fields
            WriteNameIfNotEmpty(writer, rna.Name);
            WriteRecommendedProteinNameIfNotEmpty(writer, rna.FullName);

            // Gene/organism
            WriteGeneNames(writer, rna.GeneNames);
            WriteOrganismIfNotEmpty(writer, rna.Organism);

            // Proteolysis products (no special null-begin handling here to preserve original behavior)
            WriteProteolysisProductsRna(writer, rna.TruncationProducts);

            // Base modification features
            WriteModifiedResidueFeatures(writer, rna, null, additionalMods, newModResEntries, orderModIds: false);

            // Sequence variants and their subfeatures (variant-specific mods)
            WriteRnaSequenceVariantFeatures(writer, rna, additionalMods, newModResEntries);

            // Sequence
            WriteRnaSequenceElement(writer, rna);

            writer.WriteEndElement(); // entry
        }

        /// <summary>
        /// Writes a complete protein entry with metadata, features, and sequence.
        /// </summary>
        private static void WriteProteinEntry(
            XmlWriter writer,
            Protein protein,
            bool isAppliedVariantEntry,
            bool updateTimeStamp,
            bool includeAppliedVariantFeatures,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalMods,
            Dictionary<string, int> newModResEntries)
        {
            writer.WriteStartElement("entry", "http://uniprot.org/uniprot");
            writer.WriteAttributeString("dataset", protein.DatasetEntryTag);
            writer.WriteAttributeString("created", protein.CreatedEntryTag);
            writer.WriteAttributeString("modified", updateTimeStamp ? DateTime.Now.ToString("yyyy-MM-dd") : protein.ModifiedEntryTag);
            writer.WriteAttributeString("version", protein.VersionEntryTag);

            if (isAppliedVariantEntry)
            {
                writer.WriteAttributeString("variant", "true");
            }

            // Accession and names
            WriteAccession(writer, protein.Accession);
            WriteNameIfNotNull(writer, protein.Name);
            WriteRecommendedProteinNameIfNotNull(writer, protein.FullName);

            // Gene/organism
            WriteGeneNames(writer, protein.GeneNames);
            WriteOrganismIfNotNull(writer, protein.Organism);

            // Database references
            WriteDatabaseReferences(writer, protein.DatabaseReferences);

            // Proteolysis products (with null-begin as status="unknown")
            WriteProteolysisProductsProtein(writer, protein.TruncationProducts);

            // Base modification features
            WriteModifiedResidueFeatures(writer, protein, null, additionalMods, newModResEntries, orderModIds: true);

            // Sequence variant features (candidate vs applied)
            WriteProteinSequenceVariantFeatures(writer, protein, isAppliedVariantEntry, includeAppliedVariantFeatures, additionalMods, newModResEntries);

            // Disulfide bonds
            WriteDisulfideBonds(writer, protein.DisulfideBonds);

            // Splice sites
            WriteSpliceSites(writer, protein.SpliceSites);

            // Sequence
            WriteProteinSequenceElement(writer, protein);

            writer.WriteEndElement(); // entry
        }

        /// <summary>
        /// Writes a human-readable "modified residue" feature set for a biopolymer, optionally variant-scoped.
        /// </summary>
        private static void WriteModifiedResidueFeatures(
            XmlWriter writer,
            IBioPolymer bioPolymer,
            SequenceVariation seqVar,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalMods,
            Dictionary<string, int> newModResEntries,
            bool orderModIds)
        {
            var modsForThis = GetModsForThisBioPolymer(bioPolymer, seqVar, additionalMods, newModResEntries);

            foreach (var positionModKvp in modsForThis.OrderBy(kv => kv.Key))
            {
                IEnumerable<string> ids = positionModKvp.Value;
                if (orderModIds) ids = ids.OrderBy(m => m, StringComparer.Ordinal);

                foreach (var modId in ids)
                {
                    writer.WriteStartElement("feature");
                    writer.WriteAttributeString("type", "modified residue");
                    writer.WriteAttributeString("description", modId);
                    writer.WriteStartElement("location");
                    writer.WriteStartElement(seqVar == null ? "position" : "subposition");
                    writer.WriteAttributeString(seqVar == null ? "position" : "subposition", positionModKvp.Key.ToString(CultureInfo.InvariantCulture));
                    writer.WriteEndElement(); // position/subposition
                    writer.WriteEndElement(); // location
                    writer.WriteEndElement(); // feature or subfeature
                }
            }
        }

        /// <summary>
        /// Writes RNA sequence variant features and variant-mod subfeatures. Ensures robust non-empty descriptions (VCF Description → VCF.ToString → SimpleString → synthesized).
        /// </summary>
        private static void WriteRnaSequenceVariantFeatures(
            XmlWriter writer,
            RNA rna,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalMods,
            Dictionary<string, int> newModResEntries)
        {
            foreach (var sv in (rna.SequenceVariations ?? new List<SequenceVariation>())
                .OrderBy(sv => sv.OneBasedBeginPosition)
                .ThenBy(sv => sv.VariantSequence ?? string.Empty))
            {
                if (sv == null)
                    continue;

                // Build a guaranteed non-empty description (aligned with protein logic)
                string description =
                    sv.Description ??
                    sv.VariantCallFormatData?.Description ??
                    sv.VariantCallFormatData?.ToString() ??
                    sv.SimpleString();

                if (string.IsNullOrWhiteSpace(description))
                {
                    var orig = sv.OriginalSequence ?? string.Empty;
                    var varSeq = sv.VariantSequence ?? string.Empty;
                    if (!string.IsNullOrEmpty(orig) && !string.IsNullOrEmpty(varSeq))
                    {
                        description = sv.OneBasedBeginPosition == sv.OneBasedEndPosition
                            ? $"{orig}{sv.OneBasedBeginPosition}{varSeq}"
                            : $"{orig}{sv.OneBasedBeginPosition}-{sv.OneBasedEndPosition}{varSeq}";
                    }
                    else
                    {
                        description = "sequence variant";
                    }
                }

                writer.WriteStartElement("feature");
                writer.WriteAttributeString("type", "sequence variant");
                writer.WriteAttributeString("description", description);

                writer.WriteStartElement("original");
                writer.WriteString(sv.OriginalSequence);
                writer.WriteEndElement();

                writer.WriteStartElement("variation");
                writer.WriteString(sv.VariantSequence);
                writer.WriteEndElement();

                writer.WriteStartElement("location");
                WriteSpanOrPointLocation(writer, sv.OneBasedBeginPosition, sv.OneBasedEndPosition);

                // Variant-specific modified residues as subfeatures
                foreach (var hmm in GetModsForThisBioPolymer(rna, sv, additionalMods, newModResEntries).OrderBy(b => b.Key))
                {
                    foreach (var modId in hmm.Value)
                    {
                        writer.WriteStartElement("subfeature");
                        writer.WriteAttributeString("type", "modified residue");
                        writer.WriteAttributeString("description", modId);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("subposition");
                        writer.WriteAttributeString("subposition", hmm.Key.ToString(CultureInfo.InvariantCulture));
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }
                }

                writer.WriteEndElement(); // location
                writer.WriteEndElement(); // feature
            }
        }

        /// <summary>
        /// Writes protein sequence variant features (candidate or applied) including subfeatures for variant-specific mods.
        /// Ensures a robust, non-empty description string.
        /// </summary>
        private static void WriteProteinSequenceVariantFeatures(
            XmlWriter writer,
            Protein protein,
            bool isAppliedVariantEntry,
            bool includeAppliedVariantFeatures,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalMods,
            Dictionary<string, int> newModResEntries)
        {
            IEnumerable<SequenceVariation> variantFeaturesSource =
                (protein.SequenceVariations ?? Enumerable.Empty<SequenceVariation>());

            if (isAppliedVariantEntry && includeAppliedVariantFeatures)
            {
                variantFeaturesSource = protein.AppliedSequenceVariations ?? new List<SequenceVariation>();
            }

            foreach (var sv in variantFeaturesSource
                .OrderBy(sv => sv.OneBasedBeginPosition)
                .ThenBy(sv => sv.VariantSequence ?? string.Empty))
            {
                if (sv == null) continue;

                string description =
                    sv.Description ??
                    sv.VariantCallFormatData?.Description ??
                    sv.VariantCallFormatData?.ToString() ??
                    sv.SimpleString();

                if (string.IsNullOrWhiteSpace(description))
                {
                    var orig = sv.OriginalSequence ?? string.Empty;
                    var varSeq = sv.VariantSequence ?? string.Empty;
                    if (!string.IsNullOrEmpty(orig) && !string.IsNullOrEmpty(varSeq))
                    {
                        description = sv.OneBasedBeginPosition == sv.OneBasedEndPosition
                            ? $"{orig}{sv.OneBasedBeginPosition}{varSeq}"
                            : $"{orig}{sv.OneBasedBeginPosition}-{sv.OneBasedEndPosition}{varSeq}";
                    }
                    else
                    {
                        description = "sequence variant";
                    }
                }

                writer.WriteStartElement("feature");
                writer.WriteAttributeString("type", "sequence variant");
                writer.WriteAttributeString("description", description);

                writer.WriteStartElement("original");
                writer.WriteString(sv.OriginalSequence ?? string.Empty);
                writer.WriteEndElement();

                writer.WriteStartElement("variation");
                writer.WriteString(sv.VariantSequence ?? string.Empty);
                writer.WriteEndElement();

                writer.WriteStartElement("location");
                WriteSpanOrPointLocation(writer, sv.OneBasedBeginPosition, sv.OneBasedEndPosition);

                // Variant-specific modified residues as subfeatures (ordered by mod id for stable output)
                foreach (var hmm in GetModsForThisBioPolymer(protein, sv, additionalMods, newModResEntries).OrderBy(b => b.Key))
                {
                    foreach (var modId in hmm.Value.OrderBy(m => m, StringComparer.Ordinal))
                    {
                        writer.WriteStartElement("subfeature");
                        writer.WriteAttributeString("type", "modified residue");
                        writer.WriteAttributeString("description", modId);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("subposition");
                        writer.WriteAttributeString("subposition", hmm.Key.ToString(CultureInfo.InvariantCulture));
                        writer.WriteEndElement(); // subposition
                        writer.WriteEndElement(); // location
                        writer.WriteEndElement(); // subfeature
                    }
                }

                writer.WriteEndElement(); // location
                writer.WriteEndElement(); // feature
            }
        }

        /// <summary>
        /// Writes proteolysis products for proteins; if begin is null, emits status="unknown" instead of position.
        /// </summary>
        private static void WriteProteolysisProductsProtein(XmlWriter writer, IEnumerable<TruncationProduct> products)
        {
            var proteolysisProducts = (products ?? Enumerable.Empty<TruncationProduct>())
                .Where(p => !p.Type.Contains("truncation"))
                .OrderBy(p => p.OneBasedBeginPosition)
                .ToList();

            foreach (var proteolysisProduct in proteolysisProducts)
            {
                writer.WriteStartElement("feature");
                writer.WriteAttributeString("type", proteolysisProduct.Type.Split('(')[0]);
                writer.WriteStartElement("location");
                writer.WriteStartElement("begin");

                if (proteolysisProduct.OneBasedBeginPosition == null)
                {
                    writer.WriteAttributeString("status", "unknown");
                }
                else
                {
                    writer.WriteAttributeString("position", proteolysisProduct.OneBasedBeginPosition.ToString());
                }

                writer.WriteEndElement(); // begin
                writer.WriteStartElement("end");
                writer.WriteAttributeString("position", proteolysisProduct.OneBasedEndPosition.ToString());
                writer.WriteEndElement(); // end
                writer.WriteEndElement(); // location
                writer.WriteEndElement(); // feature
            }
        }

        /// <summary>
        /// Writes proteolysis products for RNA; preserves original behavior for begin position handling.
        /// </summary>
        private static void WriteProteolysisProductsRna(XmlWriter writer, IEnumerable<TruncationProduct> products)
        {
            var proteolysisProducts = (products ?? Enumerable.Empty<TruncationProduct>())
                .Where(p => !p.Type.Contains("truncation"))
                .ToList();

            foreach (var proteolysisProduct in proteolysisProducts)
            {
                writer.WriteStartElement("feature");
                writer.WriteAttributeString("type", proteolysisProduct.Type.Split('(')[0]);
                writer.WriteStartElement("location");
                writer.WriteStartElement("begin");

                // Original RNA writer did not handle null begin specially
                writer.WriteAttributeString("position", proteolysisProduct.OneBasedBeginPosition.ToString());

                writer.WriteEndElement(); // begin
                writer.WriteStartElement("end");
                writer.WriteAttributeString("position", proteolysisProduct.OneBasedEndPosition.ToString());
                writer.WriteEndElement(); // end
                writer.WriteEndElement(); // location
                writer.WriteEndElement(); // feature
            }
        }

        /// <summary>
        /// Writes disulfide bond features with begin/end or single position.
        /// </summary>
        private static void WriteDisulfideBonds(XmlWriter writer, IEnumerable<DisulfideBond> bonds)
        {
            foreach (var bond in (bonds ?? Enumerable.Empty<DisulfideBond>()).OrderBy(b => b.OneBasedBeginPosition))
            {
                writer.WriteStartElement("feature");
                writer.WriteAttributeString("type", "disulfide bond");
                writer.WriteAttributeString("description", bond.Description);
                writer.WriteStartElement("location");
                WriteSpanOrPointLocation(writer, bond.OneBasedBeginPosition, bond.OneBasedEndPosition);
                writer.WriteEndElement(); // location
                writer.WriteEndElement(); // feature
            }
        }

        /// <summary>
        /// Writes splice site features with begin/end or single position.
        /// </summary>
        private static void WriteSpliceSites(XmlWriter writer, IEnumerable<SpliceSite> sites)
        {
            foreach (var site in (sites ?? Enumerable.Empty<SpliceSite>()).OrderBy(s => s.OneBasedBeginPosition))
            {
                writer.WriteStartElement("feature");
                writer.WriteAttributeString("type", "splice site");
                writer.WriteAttributeString("description", site.Description);
                writer.WriteStartElement("location");
                WriteSpanOrPointLocation(writer, site.OneBasedBeginPosition, site.OneBasedEndPosition);
                writer.WriteEndElement(); // location
                writer.WriteEndElement(); // feature
            }
        }

        /// <summary>
        /// Writes a span (begin/end) or a single position to the current "location" element.
        /// </summary>
        private static void WriteSpanOrPointLocation(XmlWriter writer, int begin, int end)
        {
            if (begin == end)
            {
                writer.WriteStartElement("position");
                writer.WriteAttributeString("position", begin.ToString(CultureInfo.InvariantCulture));
                writer.WriteEndElement();
            }
            else
            {
                writer.WriteStartElement("begin");
                writer.WriteAttributeString("position", begin.ToString(CultureInfo.InvariantCulture));
                writer.WriteEndElement();
                writer.WriteStartElement("end");
                writer.WriteAttributeString("position", end.ToString(CultureInfo.InvariantCulture));
                writer.WriteEndElement();
            }
        }

        /// <summary>
        /// Writes the UniProt-style sequence element with attributes for proteins.
        /// </summary>
        private static void WriteProteinSequenceElement(XmlWriter writer, Protein protein)
        {
            writer.WriteStartElement("sequence");
            writer.WriteAttributeString("length", protein.UniProtSequenceAttributes.Length.ToString(CultureInfo.InvariantCulture));
            writer.WriteAttributeString("mass", protein.UniProtSequenceAttributes.Mass.ToString(CultureInfo.InvariantCulture));
            writer.WriteAttributeString("checksum", protein.UniProtSequenceAttributes.Checksum);
            writer.WriteAttributeString("modified", protein.UniProtSequenceAttributes.EntryModified.ToString("yyyy-MM-dd"));
            writer.WriteAttributeString("version", protein.UniProtSequenceAttributes.SequenceVersion.ToString(CultureInfo.InvariantCulture));

            if (protein.UniProtSequenceAttributes.IsPrecursor != null)
            {
                writer.WriteAttributeString("precursor", protein.UniProtSequenceAttributes.IsPrecursor.Value.ToString().ToLowerInvariant());
            }

            if (protein.UniProtSequenceAttributes.Fragment != UniProtSequenceAttributes.FragmentType.unspecified)
            {
                writer.WriteAttributeString("fragment", protein.UniProtSequenceAttributes.Fragment.ToString().ToLowerInvariant());
            }

            writer.WriteString(protein.BaseSequence);
            writer.WriteEndElement(); // sequence
        }

        /// <summary>
        /// Writes the simple sequence element for RNA.
        /// </summary>
        private static void WriteRnaSequenceElement(XmlWriter writer, RNA rna)
        {
            writer.WriteStartElement("sequence");
            writer.WriteAttributeString("length", rna.Length.ToString(CultureInfo.InvariantCulture));
            writer.WriteString(rna.BaseSequence);
            writer.WriteEndElement();
        }

        /// <summary>
        /// Writes an accession element.
        /// </summary>
        private static void WriteAccession(XmlWriter writer, string accession)
        {
            writer.WriteStartElement("accession");
            writer.WriteString(accession);
            writer.WriteEndElement();
        }

        /// <summary>
        /// Writes the display name if not null.
        /// </summary>
        private static void WriteNameIfNotNull(XmlWriter writer, string name)
        {
            if (name == null) return;
            writer.WriteStartElement("name");
            writer.WriteString(name);
            writer.WriteEndElement();
        }

        /// <summary>
        /// Writes the display name if not null/empty/whitespace (RNA variant).
        /// </summary>
        private static void WriteNameIfNotEmpty(XmlWriter writer, string name)
        {
            if (!name.IsNotNullOrEmptyOrWhiteSpace()) return;
            writer.WriteStartElement("name");
            writer.WriteString(name);
            writer.WriteEndElement();
        }

        /// <summary>
        /// Writes the recommendedName/fullName block if FullName is set (protein).
        /// </summary>
        private static void WriteRecommendedProteinNameIfNotNull(XmlWriter writer, string fullName)
        {
            if (fullName == null) return;
            writer.WriteStartElement("protein");
            writer.WriteStartElement("recommendedName");
            writer.WriteStartElement("fullName");
            writer.WriteString(fullName);
            writer.WriteEndElement(); // fullName
            writer.WriteEndElement(); // recommendedName
            writer.WriteEndElement(); // protein
        }

        /// <summary>
        /// Writes the recommendedName/fullName block if FullName is not empty (RNA).
        /// </summary>
        private static void WriteRecommendedProteinNameIfNotEmpty(XmlWriter writer, string fullName)
        {
            if (!fullName.IsNotNullOrEmptyOrWhiteSpace()) return;
            writer.WriteStartElement("protein");
            writer.WriteStartElement("recommendedName");
            writer.WriteStartElement("fullName");
            writer.WriteString(fullName);
            writer.WriteEndElement(); // fullName
            writer.WriteEndElement(); // recommendedName
            writer.WriteEndElement(); // protein
        }

        /// <summary>
        /// Writes gene names.
        /// </summary>
        private static void WriteGeneNames(XmlWriter writer, IEnumerable<Tuple<string, string>> geneNames)
        {
            writer.WriteStartElement("gene");
            foreach (var geneName in (geneNames ?? Enumerable.Empty<Tuple<string, string>>()))
            {
                writer.WriteStartElement("name");
                writer.WriteAttributeString("type", geneName.Item1);
                writer.WriteString(geneName.Item2);
                writer.WriteEndElement();
            }
            writer.WriteEndElement();
        }

        /// <summary>
        /// Writes organism block if present (protein).
        /// </summary>
        private static void WriteOrganismIfNotNull(XmlWriter writer, string organism)
        {
            if (organism == null) return;
            writer.WriteStartElement("organism");
            writer.WriteStartElement("name");
            writer.WriteAttributeString("type", "scientific");
            writer.WriteString(organism);
            writer.WriteEndElement(); // name
            writer.WriteEndElement(); // organism
        }

        /// <summary>
        /// Writes organism block if string is not empty (RNA).
        /// </summary>
        private static void WriteOrganismIfNotEmpty(XmlWriter writer, string organism)
        {
            if (!organism.IsNotNullOrEmptyOrWhiteSpace()) return;
            writer.WriteStartElement("organism");
            writer.WriteStartElement("name");
            writer.WriteAttributeString("type", "scientific");
            writer.WriteString(organism);
            writer.WriteEndElement(); // name
            writer.WriteEndElement(); // organism
        }

        /// <summary>
        /// Writes database references with sorted properties for stability.
        /// </summary>
        private static void WriteDatabaseReferences(XmlWriter writer, IEnumerable<DatabaseReference> dbRefs)
        {
            foreach (var dbRef in (dbRefs ?? Enumerable.Empty<DatabaseReference>()))
            {
                writer.WriteStartElement("dbReference");
                writer.WriteAttributeString("type", dbRef.Type);
                writer.WriteAttributeString("id", dbRef.Id);

                foreach (Tuple<string, string> property in dbRef.Properties.OrderBy(t => t.Item1).ThenBy(t => t.Item2))
                {
                    writer.WriteStartElement("property");
                    writer.WriteAttributeString("type", property.Item1);
                    writer.WriteAttributeString("value", property.Item2);
                    writer.WriteEndElement();
                }

                writer.WriteEndElement();
            }
        }

        /// <summary>
        /// Returns true if a protein is an applied variant entry that should be annotated as such.
        /// </summary>
        private static bool DetermineIsAppliedVariantEntry(Protein protein, bool includeAppliedVariantEntries)
        {
            var consensus = protein.ConsensusVariant as Protein;
            return includeAppliedVariantEntries
                   && consensus != null
                   && !ReferenceEquals(protein, consensus)
                   && protein.AppliedSequenceVariations != null
                   && protein.AppliedSequenceVariations.Count > 0;
        }

        /// <summary>
        /// Creates indented XML writer settings.
        /// </summary>
        private static XmlWriterSettings CreateIndentedWriterSettings()
        {
            return new XmlWriterSettings
            {
                Indent = true,
                IndentChars = "  "
            };
        }

        /// <summary>
        /// Writes the mzLibProteinDb start element and XML declaration.
        /// </summary>
        private static void WriteStartDocument(XmlWriter writer)
        {
            writer.WriteStartDocument();
            writer.WriteStartElement("mzLibProteinDb");
        }

        /// <summary>
        /// Closes the mzLibProteinDb element and ends the document.
        /// </summary>
        private static void WriteEndDocument(XmlWriter writer)
        {
            writer.WriteEndElement(); // mzLibProteinDb
            writer.WriteEndDocument();
        }

        /// <summary>
        /// Gathers modified residue identifiers for a polymer (optionally variant-scoped), merges additional mods,
        /// and updates counts of new "modified residue" entries introduced by AdditionalMods.
        /// </summary>
        private static Dictionary<int, HashSet<string>> GetModsForThisBioPolymer(
            IBioPolymer protein,
            SequenceVariation seqvar,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins,
            Dictionary<string, int> newModResEntries)
        {
            var modsToWriteForThisSpecificProtein = new Dictionary<int, HashSet<string>>();

            // Select the appropriate modification dictionary (variant-specific if seqvar != null).
            // Each side guarantees a non-null dictionary (falls back to new Dictionary<,>()), so no further null check needed.
            var primaryModDict = seqvar == null
                ? (protein.OneBasedPossibleLocalizedModifications ?? new Dictionary<int, List<Modification>>())
                : (seqvar.OneBasedModifications ?? new Dictionary<int, List<Modification>>());

            foreach (var mods in primaryModDict)
            {
                if (mods.Value == null) continue;
                foreach (var mod in mods.Value)
                {
                    if (mod == null) continue;
                    if (modsToWriteForThisSpecificProtein.TryGetValue(mods.Key, out var set))
                        set.Add(mod.IdWithMotif);
                    else
                        modsToWriteForThisSpecificProtein.Add(mods.Key, new HashSet<string> { mod.IdWithMotif });
                }
            }

            // Additional externally supplied mods (accession changes if seqvar is applied)
            string accession = seqvar == null
                ? protein.Accession
                : VariantApplication.GetAccession(protein, new[] { seqvar });

            if (additionalModsToAddToProteins != null &&
                accession != null &&
                additionalModsToAddToProteins.TryGetValue(accession, out var extraMods))
            {
                foreach (var (pos, mod) in extraMods.Where(t => t != null))
                {
                    if (mod == null) continue;

                    bool added;
                    if (modsToWriteForThisSpecificProtein.TryGetValue(pos, out var set))
                        added = set.Add(mod.IdWithMotif);
                    else
                    {
                        modsToWriteForThisSpecificProtein.Add(pos, new HashSet<string> { mod.IdWithMotif });
                        added = true;
                    }

                    if (added)
                    {
                        if (newModResEntries.ContainsKey(mod.IdWithMotif))
                            newModResEntries[mod.IdWithMotif]++;
                        else
                            newModResEntries.Add(mod.IdWithMotif, 1);
                    }
                }
            }

            return modsToWriteForThisSpecificProtein;
        }
    }
}