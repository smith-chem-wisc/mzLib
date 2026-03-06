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
using Omics.Modifications.Conversion;
using Transcriptomics;
using System.Data;

namespace UsefulProteomicsDatabases
{

    /// <summary>
    /// Provides methods for writing rna and nucleic acid databases to XML and FASTA formats.
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
        public static Dictionary<string, int> WriteXmlDatabase(Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins, List<IBioPolymer> bioPolymerList, string outputFileName)
        {
            return WriteXmlDatabase(additionalModsToAddToProteins, bioPolymerList, outputFileName, null);
        }

        public static Dictionary<string, int> WriteXmlDatabase(
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins,
            List<IBioPolymer> bioPolymerList,
            string outputFileName,
            ProteinDbWriterConversionOptions? conversionOptions)
        {
            return bioPolymerList.Any(p => p is Protein)
                ? WriteXmlDatabase(additionalModsToAddToProteins, bioPolymerList.Cast<Protein>().ToList(), outputFileName, false, conversionOptions)
                : WriteXmlDatabase(additionalModsToAddToProteins, bioPolymerList.Cast<RNA>().ToList(), outputFileName, false, conversionOptions);
        }

        /// <summary>
        /// Writes an XML database for a list of nucleic acid sequences, including additional modifications.
        /// </summary>
        /// <param name="additionalModsToAddToNucleicAcids">A dictionary of additional modifications to add to proteins.</param>
        /// <param name="nucleicAcidList">A list of nucleic acid sequences to be written to the database.</param>
        /// <param name="outputFileName">The name of the output XML file.</param>
        /// <returns>A dictionary of new modification residue entries.</returns>
        /// <remarks>
        /// Several chunks of code are commented out. These are blocks that are intended to be implmented in the future, but
        /// are not necessary for the bare bones implementation of Transcriptomics
        /// </remarks>
        public static Dictionary<string, int> WriteXmlDatabase(Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToNucleicAcids, List<RNA> nucleicAcidList, string outputFileName, bool updateTimeStamp = false)
        {
            return WriteXmlDatabase(additionalModsToAddToNucleicAcids, nucleicAcidList, outputFileName, updateTimeStamp, null);
        }

        public static Dictionary<string, int> WriteXmlDatabase(
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToNucleicAcids,
            List<RNA> nucleicAcidList,
            string outputFileName,
            bool updateTimeStamp,
            ProteinDbWriterConversionOptions? conversionOptions)
        {
            additionalModsToAddToNucleicAcids = additionalModsToAddToNucleicAcids ?? new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var conversionContext = CreateConversionContext(conversionOptions);

            // write nonvariant rna (for cases where variants aren't applied, this just gets the rna itself)
            var nonVariantRna = nucleicAcidList.Select(p => p.ConsensusVariant).Distinct().ToList();

            var xmlWriterSettings = new XmlWriterSettings
            {
                Indent = true,
                IndentChars = "  "
            };

            Dictionary<string, int> newModResEntries = new Dictionary<string, int>();
            using (XmlWriter writer = XmlWriter.Create(outputFileName, xmlWriterSettings))
            {
                writer.WriteStartDocument();
                writer.WriteStartElement("mzLibProteinDb");

                List<Modification> myModificationList = new List<Modification>();
                foreach (var p in nonVariantRna)
                {
                    foreach (KeyValuePair<int, List<Modification>> entry in p.OneBasedPossibleLocalizedModifications)
                    {
                        myModificationList.AddRange(entry.Value);
                    }
                }

                // get modifications from nucleic acid list and concatenate the modifications discovered in GPTMDictionary
                var headerModifications = PrepareModificationDefinitions(
                    nonVariantRna
                        .SelectMany(p => p.SequenceVariations
                            .SelectMany(sv => sv.OneBasedModifications)
                            .Concat(p.OneBasedPossibleLocalizedModifications)
                            .SelectMany(kv => kv.Value))
                        .Concat(additionalModsToAddToNucleicAcids
                            .Where(kv => nonVariantRna
                                .SelectMany(p => p.SequenceVariations
                                    .Select(sv => VariantApplication.GetAccession(p, new[] { sv })).Concat(new[] { p.Accession }))
                                .Contains(kv.Key))
                            .SelectMany(kv => kv.Value.Select(v => v.Item2))),
                    conversionContext);

                foreach (Modification mod in headerModifications.OrderBy(m => m.IdWithMotif))
                {
                    writer.WriteStartElement("modification");
                    writer.WriteString(mod.ToString() + Environment.NewLine + "//");
                    writer.WriteEndElement();
                }

                foreach (var nucleicAcid in nonVariantRna)
                {
                    writer.WriteStartElement("entry", "undefined"); //this should be a website with the XSD namespace
                    //writer.WriteAttributeString("dataset", nucleicAcid.DatasetEntryTag);
                    //writer.WriteAttributeString("created", nucleicAcid.CreatedEntryTag);
                    //if (updateTimeStamp)
                    //{
                    //    writer.WriteAttributeString("modified", DateTime.Now.ToString("yyyy-MM-dd"));
                    //}
                    //else
                    //{
                    //    writer.WriteAttributeString("modified", nucleicAcid.ModifiedEntryTag);
                    //}
                    //writer.WriteAttributeString("version", nucleicAcid.VersionEntryTag);
                    writer.WriteStartElement("accession");
                    writer.WriteString(nucleicAcid.Accession);
                    writer.WriteEndElement();

                    if (nucleicAcid.Name.IsNotNullOrEmptyOrWhiteSpace())
                    {
                        writer.WriteStartElement("name");
                        writer.WriteString(nucleicAcid.Name);
                        writer.WriteEndElement();
                    }

                    if (nucleicAcid.FullName.IsNotNullOrEmptyOrWhiteSpace())
                    {
                        writer.WriteStartElement("protein");
                        writer.WriteStartElement("recommendedName");
                        writer.WriteStartElement("fullName");
                        writer.WriteString(nucleicAcid.FullName);
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    writer.WriteStartElement("gene");
                    foreach (var geneName in nucleicAcid.GeneNames)
                    {
                        writer.WriteStartElement("name");
                        writer.WriteAttributeString("type", geneName.Item1);
                        writer.WriteString(geneName.Item2);
                        writer.WriteEndElement();
                    }
                    writer.WriteEndElement();

                    if (nucleicAcid.Organism.IsNotNullOrEmptyOrWhiteSpace())
                    {
                        writer.WriteStartElement("organism");
                        writer.WriteStartElement("name");
                        writer.WriteAttributeString("type", "scientific");
                        writer.WriteString(nucleicAcid.Organism);
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    //foreach (var dbRef in nucleicAcid)
                    //{
                    //    writer.WriteStartElement("dbReference");
                    //    writer.WriteAttributeString("type", dbRef.Type);
                    //    writer.WriteAttributeString("id", dbRef.Id);
                    //    foreach (Tuple<string, string> property in dbRef.Properties)
                    //    {
                    //        writer.WriteStartElement("property");
                    //        writer.WriteAttributeString("type", property.Item1);
                    //        writer.WriteAttributeString("value", property.Item2);
                    //        writer.WriteEndElement();
                    //    }
                    //    writer.WriteEndElement();
                    //}

                    List<TruncationProduct> proteolysisProducts = nucleicAcid.TruncationProducts.Where(p => !p.Type.Contains("truncation")).ToList();
                    foreach (var proteolysisProduct in proteolysisProducts)
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", proteolysisProduct.Type.Split('(')[0]);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("begin");

                        //TODO: handle proteolysis products with null begin position
                        //see rna writer for example. 

                        writer.WriteAttributeString("position", proteolysisProduct.OneBasedBeginPosition.ToString());
                        writer.WriteEndElement();
                        writer.WriteStartElement("end");
                        writer.WriteAttributeString("position", proteolysisProduct.OneBasedEndPosition.ToString());
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    foreach (var hm in GetModsForThisBioPolymer(nucleicAcid, null, additionalModsToAddToNucleicAcids, newModResEntries, conversionContext).OrderBy(b => b.Key))
                    {
                        foreach (var modId in hm.Value)
                        {
                            writer.WriteStartElement("feature");
                            writer.WriteAttributeString("type", "modified residue");
                            writer.WriteAttributeString("description", modId);
                            writer.WriteStartElement("location");
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", hm.Key.ToString(CultureInfo.InvariantCulture));
                            writer.WriteEndElement();
                            writer.WriteEndElement();
                            writer.WriteEndElement();
                        }
                    }

                    foreach (var hm in nucleicAcid.SequenceVariations.OrderBy(sv => sv.OneBasedBeginPosition).ThenBy(sv => sv.VariantSequence))
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", "sequence variant");
                        writer.WriteAttributeString("description", hm.VariantCallFormatDataString.ToString());
                        writer.WriteStartElement("original");
                        writer.WriteString(hm.OriginalSequence);
                        writer.WriteEndElement(); // original
                        writer.WriteStartElement("variation");
                        writer.WriteString(hm.VariantSequence);
                        writer.WriteEndElement(); // variation
                        writer.WriteStartElement("location");
                        if (hm.OneBasedBeginPosition == hm.OneBasedEndPosition)
                        {
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                        }
                        else
                        {
                            writer.WriteStartElement("begin");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                            writer.WriteStartElement("end");
                            writer.WriteAttributeString("position", hm.OneBasedEndPosition.ToString());
                            writer.WriteEndElement();
                        }
                        foreach (var hmm in GetModsForThisBioPolymer(nucleicAcid, hm, additionalModsToAddToNucleicAcids, newModResEntries, conversionContext).OrderBy(b => b.Key))
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

                    //foreach (var hm in nucleicAcid.SpliceSites)
                    //{
                    //    writer.WriteStartElement("feature");
                    //    writer.WriteAttributeString("type", "splice site");
                    //    writer.WriteAttributeString("description", hm.VariantCallFormatDataString);
                    //    writer.WriteStartElement("location");
                    //    if (hm.OneBasedBeginPosition == hm.OneBasedEndPosition)
                    //    {
                    //        writer.WriteStartElement("position");
                    //        writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                    //        writer.WriteEndElement();
                    //    }
                    //    else
                    //    {
                    //        writer.WriteStartElement("begin");
                    //        writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                    //        writer.WriteEndElement();
                    //        writer.WriteStartElement("end");
                    //        writer.WriteAttributeString("position", hm.OneBasedEndPosition.ToString());
                    //        writer.WriteEndElement();
                    //    }
                    //    writer.WriteEndElement(); // location
                    //    writer.WriteEndElement(); // feature
                    //}

                    writer.WriteStartElement("sequence");
                    writer.WriteAttributeString("length", nucleicAcid.Length.ToString(CultureInfo.InvariantCulture));
                    writer.WriteString(nucleicAcid.BaseSequence);
                    writer.WriteEndElement(); // sequence
                    writer.WriteEndElement(); // entry
                }

                writer.WriteEndElement(); // mzLibProteinDb
                writer.WriteEndDocument();
            }
            return newModResEntries;
        }

        /// <summary>
        /// Writes a rna database in mzLibProteinDb format, with additional modifications from the AdditionalModsToAddToProteins list.
        /// </summary>
        /// <param name="additionalModsToAddToProteins"></param>
        /// <param name="proteinList"></param>
        /// <param name="outputFileName"></param>
        /// <returns>The new "modified residue" entries that are added due to being in the Mods dictionary</returns>
        public static Dictionary<string, int> WriteXmlDatabase(Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins, List<Protein> proteinList, string outputFileName, bool updateTimeStamp = false)
        {
            return WriteXmlDatabase(additionalModsToAddToProteins, proteinList, outputFileName, updateTimeStamp, null);
        }

        public static Dictionary<string, int> WriteXmlDatabase(
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins,
            List<Protein> proteinList,
            string outputFileName,
            bool updateTimeStamp,
            ProteinDbWriterConversionOptions? conversionOptions)
        {
            additionalModsToAddToProteins = additionalModsToAddToProteins ?? new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var conversionContext = CreateConversionContext(conversionOptions);

            // write nonvariant proteins (for cases where variants aren't applied, this just gets the rna itself)
            var nonVariantProteins = proteinList.Select(p => p.ConsensusVariant).Distinct().ToList();

            var xmlWriterSettings = new XmlWriterSettings
            {
                Indent = true,
                IndentChars = "  "
            };

            Dictionary<string, int> newModResEntries = new Dictionary<string, int>();

            using (XmlWriter writer = XmlWriter.Create(outputFileName, xmlWriterSettings))
            {
                writer.WriteStartDocument();
                writer.WriteStartElement("mzLibProteinDb");

                List<Modification> myModificationList = new List<Modification>();
                foreach (Protein p in nonVariantProteins)
                {
                    foreach (KeyValuePair<int, List<Modification>> entry in p.OneBasedPossibleLocalizedModifications)
                    {
                        myModificationList.AddRange(entry.Value);
                    }
                }

                var headerModifications = PrepareModificationDefinitions(
                    nonVariantProteins
                        .SelectMany(p => p.SequenceVariations
                            .SelectMany(sv => sv.OneBasedModifications)
                            .Concat(p.OneBasedPossibleLocalizedModifications)
                            .SelectMany(kv => kv.Value))
                        .Concat(additionalModsToAddToProteins
                            .Where(kv => nonVariantProteins
                                .SelectMany(p => p.SequenceVariations
                                    .Select(sv => VariantApplication.GetAccession(p, new[] { sv })).Concat(new[] { p.Accession }))
                                .Contains(kv.Key))
                            .SelectMany(kv => kv.Value.Select(v => v.Item2))),
                    conversionContext);

                foreach (Modification mod in headerModifications.OrderBy(m => m.IdWithMotif))
                {
                    writer.WriteStartElement("modification");
                    writer.WriteString(mod.ToString() + Environment.NewLine + "//");
                    writer.WriteEndElement();
                }

                foreach (Protein protein in nonVariantProteins)
                {
                    writer.WriteStartElement("entry", "http://uniprot.org/uniprot");
                    writer.WriteAttributeString("dataset", protein.DatasetEntryTag);
                    writer.WriteAttributeString("created", protein.CreatedEntryTag);
                    if (updateTimeStamp)
                    {
                        writer.WriteAttributeString("modified", DateTime.Now.ToString("yyyy-MM-dd"));
                    }
                    else
                    {
                        writer.WriteAttributeString("modified", protein.ModifiedEntryTag);
                    }         
                    writer.WriteAttributeString("version", protein.VersionEntryTag);
                    writer.WriteStartElement("accession");
                    writer.WriteString(protein.Accession);
                    writer.WriteEndElement();

                    if (protein.Name != null)
                    {
                        writer.WriteStartElement("name");
                        writer.WriteString(protein.Name);
                        writer.WriteEndElement();
                    }

                    if (protein.FullName != null)
                    {
                        writer.WriteStartElement("protein");
                        writer.WriteStartElement("recommendedName");
                        writer.WriteStartElement("fullName");
                        writer.WriteString(protein.FullName);
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    writer.WriteStartElement("gene");
                    foreach (var gene_name in protein.GeneNames)
                    {
                        writer.WriteStartElement("name");
                        writer.WriteAttributeString("type", gene_name.Item1);
                        writer.WriteString(gene_name.Item2);
                        writer.WriteEndElement();
                    }
                    writer.WriteEndElement();

                    if (protein.Organism != null)
                    {
                        writer.WriteStartElement("organism");
                        writer.WriteStartElement("name");
                        writer.WriteAttributeString("type", "scientific");
                        writer.WriteString(protein.Organism);
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    foreach (var dbRef in protein.DatabaseReferences)
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

                    //for now we are not going to write top-down truncations generated for top-down truncation search. 
                    //some day we could write those if observed
                    //the truncation designation is contained in the "type" field of TruncationProduct
                    List<TruncationProduct> proteolysisProducts = protein.TruncationProducts.Where(p => !p.Type.Contains("truncation"))
                        .OrderBy(p => p.OneBasedBeginPosition).ToList();
                    foreach (var proteolysisProduct in proteolysisProducts)
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", proteolysisProduct.Type.Split('(')[0]);
                        writer.WriteStartElement("location");
                        writer.WriteStartElement("begin");

                        if(proteolysisProduct.OneBasedBeginPosition == null)
                        {
                            writer.WriteAttributeString("status", "unknown");
                        }
                        else
                        {
                            writer.WriteAttributeString("position", proteolysisProduct.OneBasedBeginPosition.ToString());
                        }

                        //writer.WriteAttributeString("position", proteolysisProduct.OneBasedBeginPosition.ToString());
                        writer.WriteEndElement();
                        writer.WriteStartElement("end");
                        writer.WriteAttributeString("position", proteolysisProduct.OneBasedEndPosition.ToString());
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                        writer.WriteEndElement();
                    }

                    foreach (var positionModKvp in GetModsForThisBioPolymer(protein, null, additionalModsToAddToProteins, newModResEntries, conversionContext).OrderBy(b => b.Key))
                    {
                        foreach (var modId in positionModKvp.Value.OrderBy(mod => mod))
                        {
                            writer.WriteStartElement("feature");
                            writer.WriteAttributeString("type", "modified residue");
                            writer.WriteAttributeString("description", modId);
                            writer.WriteStartElement("location");
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", positionModKvp.Key.ToString(CultureInfo.InvariantCulture));
                            writer.WriteEndElement();
                            writer.WriteEndElement();
                            writer.WriteEndElement();
                        }
                    }

                    
                    foreach (var hm in protein.SequenceVariations.OrderBy(sv => sv.OneBasedBeginPosition).ThenBy(sv => sv.VariantSequence)) 
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", "sequence variant");
                        writer.WriteAttributeString("description", hm.VariantCallFormatDataString.ToString());
                        writer.WriteStartElement("original");
                        writer.WriteString(hm.OriginalSequence);
                        writer.WriteEndElement(); // original
                        writer.WriteStartElement("variation");
                        writer.WriteString(hm.VariantSequence);
                        writer.WriteEndElement(); // variation
                        writer.WriteStartElement("location");
                        if (hm.OneBasedBeginPosition == hm.OneBasedEndPosition)
                        {
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                        }
                        else
                        {
                            writer.WriteStartElement("begin");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                            writer.WriteStartElement("end");
                            writer.WriteAttributeString("position", hm.OneBasedEndPosition.ToString());
                            writer.WriteEndElement();
                        }
                        foreach (var hmm in GetModsForThisBioPolymer(protein, hm, additionalModsToAddToProteins, newModResEntries, conversionContext).OrderBy(b => b.Key))
                        {
                            foreach (var modId in hmm.Value.OrderBy(mod => mod))
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

                    foreach (var hm in protein.DisulfideBonds.OrderBy(bond => bond.OneBasedBeginPosition))
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", "disulfide bond");
                        writer.WriteAttributeString("description", hm.Description);
                        writer.WriteStartElement("location");
                        if (hm.OneBasedBeginPosition == hm.OneBasedEndPosition)
                        {
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                        }
                        else
                        {
                            writer.WriteStartElement("begin");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                            writer.WriteStartElement("end");
                            writer.WriteAttributeString("position", hm.OneBasedEndPosition.ToString());
                            writer.WriteEndElement();
                        }
                        writer.WriteEndElement(); // location
                        writer.WriteEndElement(); // feature
                    }

                    foreach (var hm in protein.SpliceSites.OrderBy(site => site.OneBasedBeginPosition))
                    {
                        writer.WriteStartElement("feature");
                        writer.WriteAttributeString("type", "splice site");
                        writer.WriteAttributeString("description", hm.Description);
                        writer.WriteStartElement("location");
                        if (hm.OneBasedBeginPosition == hm.OneBasedEndPosition)
                        {
                            writer.WriteStartElement("position");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                        }
                        else
                        {
                            writer.WriteStartElement("begin");
                            writer.WriteAttributeString("position", hm.OneBasedBeginPosition.ToString());
                            writer.WriteEndElement();
                            writer.WriteStartElement("end");
                            writer.WriteAttributeString("position", hm.OneBasedEndPosition.ToString());
                            writer.WriteEndElement();
                        }
                        writer.WriteEndElement(); // location
                        writer.WriteEndElement(); // feature
                    }

                    writer.WriteStartElement("sequence");
                    writer.WriteAttributeString("length", protein.UniProtSequenceAttributes.Length.ToString(CultureInfo.InvariantCulture));
                    writer.WriteAttributeString("mass", protein.UniProtSequenceAttributes.Mass.ToString(CultureInfo.InvariantCulture));
                    writer.WriteAttributeString("checksum", protein.UniProtSequenceAttributes.Checksum);
                    writer.WriteAttributeString("modified", protein.UniProtSequenceAttributes.EntryModified.ToString("yyyy-MM-dd"));
                    writer.WriteAttributeString("version", protein.UniProtSequenceAttributes.SequenceVersion.ToString(CultureInfo.InvariantCulture));
                    //optional attributes
                    if (protein.UniProtSequenceAttributes.IsPrecursor != null)
                    {
                        writer.WriteAttributeString("precursor", protein.UniProtSequenceAttributes.IsPrecursor.Value.ToString().ToLowerInvariant());
                    }
                    if(protein.UniProtSequenceAttributes.Fragment != UniProtSequenceAttributes.FragmentType.unspecified)
                    {
                        writer.WriteAttributeString("fragment", protein.UniProtSequenceAttributes.Fragment.ToString().ToLowerInvariant());
                    }
                    //end optional attributes
                    writer.WriteString(protein.BaseSequence);
                    writer.WriteEndElement(); // sequence
                    writer.WriteEndElement(); // entry
                }

                writer.WriteEndElement(); // mzLibProteinDb
                writer.WriteEndDocument();
            }
            return newModResEntries;
        }


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

        public static void WriteFastaDatabase(List<RNA> rnaList, string outputFileName)
        {
            using (StreamWriter writer = new StreamWriter(outputFileName))
            {
                foreach (RNA rna in rnaList)
                {
                    var n = rna.GeneNames.FirstOrDefault();
                    string geneName = n == null ? "" : n.Item2;


                    writer.WriteLine(">mz|{0}|{1} {2} OS={3} GN={4}", rna.Accession, rna.Name, rna.FullName, rna.Organism, geneName);
                    writer.WriteLine(rna.BaseSequence);
                }
            }
        }

        private static Dictionary<int, HashSet<string>> GetModsForThisBioPolymer(
            IBioPolymer protein,
            SequenceVariation seqvar,
            Dictionary<string, HashSet<Tuple<int, Modification>>> additionalModsToAddToProteins,
            Dictionary<string, int> newModResEntries,
            SequenceConversionContext? conversionContext)
        {
            var modsToWriteForThisSpecificProtein = new Dictionary<int, HashSet<string>>();

            var primaryModDict = seqvar == null ? protein.OneBasedPossibleLocalizedModifications : seqvar.OneBasedModifications;
            var baseSequence = seqvar?.VariantSequence ?? protein.BaseSequence;
            foreach (var mods in primaryModDict)
            {
                foreach (var mod in mods.Value)
                {
                    if (!TryConvertModificationForWriter(baseSequence, mods.Key, mod, conversionContext, out var descriptor))
                    {
                        continue;
                    }

                    if (modsToWriteForThisSpecificProtein.TryGetValue(mods.Key, out HashSet<string> val))
                        val.Add(descriptor);
                    else
                        modsToWriteForThisSpecificProtein.Add(mods.Key, new HashSet<string> { descriptor });
                }
            }

            string accession = seqvar == null ? protein.Accession : VariantApplication.GetAccession(protein, new[] { seqvar }); 
            if (additionalModsToAddToProteins.ContainsKey(accession))
            {
                foreach (var ye in additionalModsToAddToProteins[accession])
                {
                    int additionalModResidueIndex = ye.Item1;
                    var additionalMod = ye.Item2;
                    if (!TryConvertModificationForWriter(baseSequence, additionalModResidueIndex, additionalMod, conversionContext, out var descriptor))
                    {
                        continue;
                    }

                    bool modAdded = false;

                    // If we already have modifications that need to be written to the specific residue, get the hash set of those mods
                    if (modsToWriteForThisSpecificProtein.TryGetValue(additionalModResidueIndex, out HashSet<string> val))
                    {
                        // Try to add the new mod to that hash set. If it's not there, modAdded=true, and it is added.
                        modAdded = val.Add(descriptor);
                    }

                    // Otherwise, no modifications currently need to be written to the residue at residueIndex, so need to create new hash set for that residue
                    else
                    {
                        modsToWriteForThisSpecificProtein.Add(additionalModResidueIndex, new HashSet<string> { descriptor });
                        modAdded = true;
                    }

                    // Finally, if a new modification has in fact been deemed worthy of being added to the database, mark that in the output dictionary
                    if (modAdded)
                    {
                        if (newModResEntries.ContainsKey(descriptor))
                        {
                            newModResEntries[descriptor]++;
                        }
                        else
                        {
                            newModResEntries.Add(descriptor, 1);
                        }
                    }
                }
            }
            return modsToWriteForThisSpecificProtein;
        }

        private static IEnumerable<Modification> PrepareModificationDefinitions(
            IEnumerable<Modification> modifications,
            SequenceConversionContext? conversionContext)
        {
            if (conversionContext == null)
            {
                return modifications.Distinct();
            }

            var converter = conversionContext.Converter;
            var prepared = new List<Modification>();

            foreach (var mod in modifications)
            {
                try
                {
                    if (converter.TryConvertModificationDefinition(mod, conversionContext.TargetConvention, out var converted) && converted != null)
                    {
                        prepared.Add(converted);
                    }
                    else if (conversionContext.HandlingMode == SequenceConversionHandlingMode.RemoveIncompatibleMods)
                    {
                        continue;
                    }
                    else if (conversionContext.HandlingMode == SequenceConversionHandlingMode.KeepOriginalAnnotation)
                    {
                        prepared.Add(mod);
                    }
                    else
                    {
                        throw new SequenceConversionException(
                            $"Unable to convert modification {mod.IdWithMotif} to {conversionContext.TargetConvention}.",
                            SequenceConversionFailureReason.NoEquivalent,
                            mod,
                            conversionContext.TargetConvention);
                    }
                }
                catch (SequenceConversionException)
                {
                    if (conversionContext.HandlingMode == SequenceConversionHandlingMode.RemoveIncompatibleMods)
                    {
                        continue;
                    }

                    if (conversionContext.HandlingMode == SequenceConversionHandlingMode.KeepOriginalAnnotation)
                    {
                        prepared.Add(mod);
                        continue;
                    }

                    throw;
                }
            }

            return prepared.Distinct();
        }

        private static bool TryConvertModificationForWriter(
            string baseSequence,
            int residueOneBasedIndex,
            Modification modification,
            SequenceConversionContext? conversionContext,
            out string descriptor)
        {
            if (conversionContext == null)
            {
                descriptor = modification.IdWithMotif;
                return true;
            }

            var converter = conversionContext.Converter;
            var key = GetSequenceConverterKey(residueOneBasedIndex, baseSequence.Length, modification);
            var success = converter.TryConvertModification(
                modification,
                baseSequence,
                key,
                conversionContext.TargetConvention,
                conversionContext.HandlingMode,
                out var converted,
                out var reason);

            if (!success)
            {
                throw new SequenceConversionException(
                    $"Unable to convert modification {modification.IdWithMotif} at position {residueOneBasedIndex} on sequence {baseSequence}.",
                    reason ?? SequenceConversionFailureReason.NoEquivalent,
                    modification,
                    conversionContext.TargetConvention);
            }

            if (converted == null)
            {
                descriptor = string.Empty;
                return false;
            }

            descriptor = converted.IdWithMotif;
            return true;
        }

        private static int GetSequenceConverterKey(int residueOneBasedIndex, int sequenceLength, Modification modification)
        {
            if (IsNTerminalRestriction(modification.LocationRestriction))
            {
                return 1;
            }

            if (IsCTerminalRestriction(modification.LocationRestriction))
            {
                return sequenceLength + 2;
            }

            return residueOneBasedIndex + 1;
        }

        private static bool IsNTerminalRestriction(string restriction)
        {
            if (string.IsNullOrEmpty(restriction))
            {
                return false;
            }

            return restriction.Contains("N-terminal", StringComparison.OrdinalIgnoreCase)
                || restriction.Contains("5'-terminal", StringComparison.OrdinalIgnoreCase)
                || restriction.Contains("Peptide N-terminal", StringComparison.OrdinalIgnoreCase)
                || restriction.Contains("Oligo 5'-terminal", StringComparison.OrdinalIgnoreCase);
        }

        private static bool IsCTerminalRestriction(string restriction)
        {
            if (string.IsNullOrEmpty(restriction))
            {
                return false;
            }

            return restriction.Contains("C-terminal", StringComparison.OrdinalIgnoreCase)
                || restriction.Contains("3'-terminal", StringComparison.OrdinalIgnoreCase)
                || restriction.Contains("Peptide C-terminal", StringComparison.OrdinalIgnoreCase)
                || restriction.Contains("Oligo 3'-terminal", StringComparison.OrdinalIgnoreCase);
        }

        private static SequenceConversionContext? CreateConversionContext(ProteinDbWriterConversionOptions? options)
        {
            if (options is not { Enabled: true })
            {
                return null;
            }

            if (options.HandlingMode is SequenceConversionHandlingMode.UsePrimarySequence
                or SequenceConversionHandlingMode.UseMassShifts
                or SequenceConversionHandlingMode.ReturnNull)
            {
                throw new NotSupportedException("ProteinDbWriter conversion only supports RemoveIncompatibleMods, KeepOriginalAnnotation, or ThrowException handling modes.");
            }

            return new SequenceConversionContext(
                options.Converter ?? SequenceConverter.Default,
                options.TargetConvention,
                options.HandlingMode);
        }

        private sealed class SequenceConversionContext
        {
            public SequenceConversionContext(
                ISequenceConverter converter,
                ModificationNamingConvention targetConvention,
                SequenceConversionHandlingMode handlingMode)
            {
                Converter = converter;
                TargetConvention = targetConvention;
                HandlingMode = handlingMode;
            }

            public ISequenceConverter Converter { get; }
            public ModificationNamingConvention TargetConvention { get; }
            public SequenceConversionHandlingMode HandlingMode { get; }
        }
    }
}
