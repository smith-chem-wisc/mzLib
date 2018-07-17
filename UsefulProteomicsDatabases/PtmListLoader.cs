using Chemistry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    public static class PtmListLoader
    {
        private static readonly Dictionary<string, char> aminoAcidCodes;

        static PtmListLoader()
        {
            aminoAcidCodes = new Dictionary<string, char>
            {
                { "Alanine", 'A' },
                { "Arginine", 'R' },
                { "Asparagine", 'N' },
                { "Aspartate", 'D' },
                { "Aspartic Acid", 'D' },
                { "Cysteine", 'C' },
                { "Glutamate", 'E' },
                { "Glutamic Acid", 'E' },
                { "Glutamine", 'Q' },
                { "Glycine", 'G' },
                { "Histidine", 'H' },
                { "Isoleucine", 'I' },
                { "Leucine", 'L' },
                { "Lysine", 'K' },
                { "Methionine", 'M' },
                { "Phenylalanine", 'F' },
                { "Proline", 'P' },
                { "Serine", 'S' },
                { "Threonine", 'T' },
                { "Tryptophan", 'W' },
                { "Tyrosine", 'Y' },
                { "Valine", 'V' }
            };
        }

        public static IEnumerable<Modification> ReadModsFromFile(string ptmListLocation)
        {
            return ReadModsFromFile(ptmListLocation, new Dictionary<string, int>()).OrderBy(b => b.id);
        }

        /// <summary>
        /// Reads a list of modifications from a text file.
        /// </summary>
        /// <param name="ptmListLocation"></param>
        /// <returns></returns>
        public static IEnumerable<Modification> ReadModsFromFile(string ptmListLocation, Dictionary<string, int> formalChargesDictionary)
        {
            using (StreamReader uniprot_mods = new StreamReader(ptmListLocation))
            {
                List<string> modification_specification = new List<string>();

                //This block will read one complete modification entry at a time until the EOF is reached.
                while (uniprot_mods.Peek() != -1) //The Peek method returns an integer value in order to determine whether the end of the file, or another error has occurred.
                {
                    string line = uniprot_mods.ReadLine();
                    modification_specification.Add(line);
                    if (line.StartsWith("//"))
                    {
                        foreach (var mod in ReadMod(modification_specification, formalChargesDictionary))
                        { yield return mod; }
                        modification_specification = new List<string>();
                    }
                }
            }
        }

        /// <summary>
        /// Reads a list of modifications from a string representation of a ptmlist text file.
        /// </summary>
        /// <param name="storedModifications"></param>
        /// <returns></returns>
        public static IEnumerable<Modification> ReadModsFromString(string storedModifications)
        {
            using (StringReader uniprot_mods = new StringReader(storedModifications))
            {
                List<string> modification_specification = new List<string>();

                while (uniprot_mods.Peek() != -1)
                {
                    string line = uniprot_mods.ReadLine();
                    modification_specification.Add(line);
                    if (line.StartsWith("//"))
                    {
                        foreach (var mod in ReadMod(modification_specification, new Dictionary<string, int>()))
                            yield return mod;
                        modification_specification = new List<string>();
                    }
                }
            }
        }

        /// <summary>
        /// Get a ModificationWithLocation from string representations of a modification specification. Returns null if the string representation is not recognized.
        /// </summary>
        /// <param name="specification"></param>
        /// <returns></returns>
        private static IEnumerable<Modification> ReadMod(List<string> specification, Dictionary<string, int> formalChargesDictionary)
        {
            // UniProt-specific fields
            string uniprotAC = null;
            string uniprotFT = null;

            // Other fields
            string id = null;
            List<string> motifs = null;
            string terminusLocalizationString = null;
            ChemicalFormula correctionFormula = null;
            double? monoisotopicMass = null;
            var externalDatabaseLinks = new Dictionary<string, IList<string>>();
            List<string> keywords = null;

            // Custom fields
            List<double> neutralLosses = null;
            List<double> diagnosticIons = null;
            string modificationType = null;

            foreach (string line in specification)
            {
                if (line.Length >= 2)
                {
                    string modKey = line.Substring(0, 2);
                    string modValue = null;
                    if (line.Length > 5)
                    {
                        try
                        {
                            modValue = line.Split('#')[0].Trim().Substring(5);
                        }
                        catch
                        {
                            //This catches a bug where there is a correct two letter code entry but no information that follows. so, when trim get's at it, the string is not at least 5 characters and then there is a crash.
                        }
                    }

                    switch (modKey)
                    {
                        case "ID": // Mandatory
                            id = modValue;
                            break;

                        case "AC": // Do not use! Only present in UniProt ptmlist
                            uniprotAC = modValue;
                            break;

                        case "FT": // Optional
                            uniprotFT = modValue;
                            break;

                        case "TG": // Which amino acid(s) or motifs is the modification on
                            motifs = new List<string>(modValue.TrimEnd('.').Split(new string[] { " or " }, StringSplitOptions.None));
                            break;

                        case "PP": // Terminus localization
                            terminusLocalizationString = modValue;
                            break;

                        case "CF": // Correction formula
                            correctionFormula = ChemicalFormula.ParseFormula(modValue.Replace(" ", string.Empty));
                            break;

                        case "MM": // Monoisotopic mass difference. Might not precisely correspond to formula!
                            {
                                if (!double.TryParse(modValue, NumberStyles.Any, CultureInfo.InvariantCulture, out double thisMM))
                                { throw new MzLibException(modValue + " is not a valid monoisotopic mass"); }
                                monoisotopicMass = thisMM;
                            }
                            break;

                        case "DR": // External database links!
                            {
                                var splitString = modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);
                                if (externalDatabaseLinks.TryGetValue(splitString[0], out IList<string> val))
                                    val.Add(splitString[1]);
                                else
                                    externalDatabaseLinks.Add(splitString[0], new List<string> { splitString[1] });
                            }
                            break;

                        case "KW": // ; Separated keywords
                            {
                                keywords = new List<string>(modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None));
                            }
                            break;

                        // NOW CUSTOM FIELDS:

                        case "NL": // Netural Losses. If field doesn't exist, single equal to 0
                            try
                            {
                                neutralLosses = new List<double>(modValue.Split(new string[] { " or " }, StringSplitOptions.RemoveEmptyEntries).Select(b => ChemicalFormula.ParseFormula(b).MonoisotopicMass));
                            }
                            catch (MzLibException)
                            {
                                neutralLosses = new List<double>(modValue.Split(new string[] { " or " }, StringSplitOptions.RemoveEmptyEntries).Select(b => double.Parse(b, CultureInfo.InvariantCulture)));
                            }
                            break;

                        case "DI": // Masses of diagnostic ions. Might just be "DI"!!! If field doesn't exist, create an empty list!
                            try
                            {
                                var ok = line;
                                var ok2 = ok.Split('#');
                                var ok3 = ok2[0];
                                var ok4 = ok3.Trim();
                                var ok5 = ok4.Substring(5);

                                diagnosticIons = new List<double>(modValue.Split(new string[] { " or " }, StringSplitOptions.RemoveEmptyEntries).Select(b => ChemicalFormula.ParseFormula(b).MonoisotopicMass));
                            }
                            catch (MzLibException)
                            {
                                diagnosticIons = new List<double>(modValue.Split(new string[] { " or " }, StringSplitOptions.RemoveEmptyEntries).Select(b => double.Parse(b, CultureInfo.InvariantCulture)));
                            }
                            break;

                        case "MT": // Modification Type. If the field doesn't exist, set to the database name
                            modificationType = modValue;
                            break;

                        case "//":
                            if (id == null)
                                throw new MzLibException("id is null");
                            if ("CROSSLNK".Equals(uniprotFT)) // Ignore crosslinks
                                break;
                            if (uniprotAC != null)
                            {
                                modificationType = "UniProt";
                                externalDatabaseLinks.Add("UniProt", new List<string> { uniprotAC });
                            }
                            if (modificationType == null)
                                throw new MzLibException("modificationType of " + id + " is null");
                            if (!monoisotopicMass.HasValue && correctionFormula != null)
                                monoisotopicMass = correctionFormula.MonoisotopicMass;

                            foreach (var dbAndAccession in externalDatabaseLinks.SelectMany(b => b.Value.Select(c => b.Key + "; " + c)))
                                if (formalChargesDictionary.ContainsKey(dbAndAccession))
                                {
                                    if (monoisotopicMass.HasValue)
                                        monoisotopicMass -= formalChargesDictionary[dbAndAccession] * Constants.ProtonMass;
                                    if (correctionFormula != null)
                                        correctionFormula.Remove(PeriodicTable.GetElement("H"), formalChargesDictionary[dbAndAccession]);
                                    break;
                                }
                            if (terminusLocalizationString == null || motifs == null)
                                yield return new Modification(id, modificationType);
                            else if (ModificationWithLocation.terminusLocalizationTypeCodes.TryGetValue(terminusLocalizationString, out TerminusLocalization terminusLocalization))
                            {
                                foreach (var singleTarget in motifs)
                                {
                                    string theMotif;
                                    if (aminoAcidCodes.TryGetValue(singleTarget, out char possibleMotifChar))
                                        theMotif = possibleMotifChar.ToString();
                                    else
                                        theMotif = singleTarget;
                                    if (ModificationMotif.TryGetMotif(theMotif, out ModificationMotif motif))
                                    {
                                        var idToUse = id;
                                        // Augment id if mulitple motifs!
                                        // Add id to keywords
                                        if (motifs.Count != 1)
                                        {
                                            if (keywords == null)
                                                keywords = new List<string> { id };
                                            else
                                                keywords.Add(id);
                                            idToUse += " on " + motif;
                                        }

                                        // Add the modification!

                                        if (!monoisotopicMass.HasValue)
                                        {
                                            // Return modification
                                            yield return new ModificationWithLocation(idToUse, modificationType, motif, terminusLocalization, externalDatabaseLinks, keywords);
                                        }
                                        else
                                        {
                                            if (correctionFormula == null)
                                            {
                                                // Return modification with mass
                                                yield return new ModificationWithMass(idToUse, modificationType, motif, terminusLocalization, monoisotopicMass.Value, externalDatabaseLinks,
                                                    keywords,
                                                    neutralLosses,
                                                    diagnosticIons);
                                            }
                                            else
                                            {
                                                // Return modification with complete information!
                                                yield return new ModificationWithMassAndCf(idToUse, modificationType, motif, terminusLocalization, correctionFormula, monoisotopicMass.Value, externalDatabaseLinks, keywords,
                                                    neutralLosses,
                                                    diagnosticIons);
                                            }
                                        }
                                    }
                                    else
                                        throw new MzLibException("Could not get motif from " + singleTarget);
                                }
                            }
                            else
                                throw new MzLibException("Could not get modification site from " + terminusLocalizationString);
                            break;
                        default:
                            break; 
                    }
                }
            }
        }
    }
}