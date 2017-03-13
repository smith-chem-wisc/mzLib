using Chemistry;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    public static class PtmListLoader
    {

        #region Private Fields

        private static readonly Dictionary<string, char> aminoAcidCodes;

        #endregion Private Fields

        #region Public Constructors

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

        #endregion Public Constructors

        #region Public Methods

        /// <summary>
        /// Reads a list of modifications from a text file.
        /// </summary>
        /// <param name="ptmListLocation"></param>
        /// <returns></returns>
        public static IEnumerable<ModificationWithLocation> ReadModsFromFile(string ptmListLocation)
        {
            using (StreamReader uniprot_mods = new StreamReader(ptmListLocation))
            {
                List<string> modification_specification = new List<string>();

                while (uniprot_mods.Peek() != -1)
                {
                    string line = uniprot_mods.ReadLine();
                    if (line.Length >= 2)
                    {
                        switch (line.Substring(0, 2))
                        {
                            case "ID":
                            case "AC":
                            case "FT": // MOD_RES CROSSLNK LIPID
                            case "TG": // Which amino acid(s) or motifs is the modification on
                            case "PP": // Terminus localization
                            case "CF": // Correction formula
                            case "MM": // Monoisotopic mass difference. Might not precisely correspond to formula!
                            case "DR": // External database links!

                            // NOW CUSTOM FIELDS:
                            case "NL": // Netural Losses. If field doesn't exist, single equal to 0
                            case "OM": // What masses are seen in histogram. If field doesn't exist, single equal to MM
                            case "DI": // Masses of diagnostic ions
                            case "MT": // Modification Type. If the field doesn't exist, set to the database name
                                modification_specification.Add(line);
                                break;

                            case "//":
                                modification_specification.Add(line);
                                ModificationWithLocation m = ReadMod(modification_specification);
                                if (m != null) yield return m;
                                modification_specification = new List<string>();
                                break;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Reads a list of modifications from a string representation of a ptmlist text file.
        /// </summary>
        /// <param name="storedModifications"></param>
        /// <returns></returns>
        public static IEnumerable<ModificationWithLocation> ReadModsFromString(string storedModifications)
        {
            using (StringReader uniprot_mods = new StringReader(storedModifications))
            {
                List<string> modification_specification = new List<string>();

                while (uniprot_mods.Peek() != -1)
                {
                    string line = uniprot_mods.ReadLine();
                    if (line.Length >= 2)
                    {
                        switch (line.Substring(0, 2))
                        {
                            case "ID":
                            case "AC":
                            case "FT": // MOD_RES CROSSLNK LIPID
                            case "TG": // Which amino acid(s) or motifs is the modification on
                            case "PP": // Terminus localization
                            case "CF": // Correction formula
                            case "MM": // Monoisotopic mass difference. Might not precisely correspond to formula!
                            case "DR": // External database links!

                            // NOW CUSTOM FIELDS:
                            case "NL": // Netural Losses. If field doesn't exist, single equal to 0
                            case "OM": // What masses are seen in histogram. If field doesn't exist, single equal to MM
                            case "DI": // Masses of diagnostic ions
                            case "MT": // Modification Type. If the field doesn't exist, set to the database name
                                modification_specification.Add(line);
                                break;

                            case "//":
                                modification_specification.Add(line);
                                ModificationWithLocation m = ReadMod(modification_specification);
                                if (m != null) yield return m;
                                modification_specification = new List<string>();
                                break;
                        }
                    }
                }
            }
        }

        #endregion Public Methods

        #region Private Methods

        /// <summary>
        /// Get a ModificationWithLocation from string representations of a modification specification. Returns null if the string representation is not recognized.
        /// </summary>
        /// <param name="specification"></param>
        /// <returns></returns>
        private static ModificationWithLocation ReadMod(List<string> specification)
        {
            // UniProt fields
            string uniprotID = null;
            Tuple<string, string> uniprotAC = null;
            string uniprotFT = null;
            IEnumerable<string> uniprotTG = null;
            string uniprotPP = null;
            ChemicalFormula uniprotCF = null;
            double? uniprotMM = null;
            var uniprotDR = new Dictionary<string, IList<string>>();

            // Custom fields
            IEnumerable<double> neutralLosses = null;
            IEnumerable<double> massesObserved = null;
            IEnumerable<double> diagnosticIons = null;
            string modificationType = null;

            ModificationWithLocation result = null;

            foreach (string line in specification)
            {
                if (line.Length >= 2)
                {
                    switch (line.Substring(0, 2))
                    {
                        case "ID":
                            uniprotID = line.Substring(5);
                            break;

                        case "AC":
                            uniprotAC = new Tuple<string, string>("uniprot", line.Substring(5));
                            break;

                        case "FT": // MOD_RES CROSSLNK LIPID
                            uniprotFT = line.Substring(5);
                            break;

                        case "TG": // Which amino acid(s) or motifs is the modification on
                            uniprotTG = new List<string>(line.Substring(5).TrimEnd('.').Split(new string[] { " or " }, StringSplitOptions.None));
                            break;

                        case "PP": // Terminus localization
                            uniprotPP = line.Substring(5);
                            break;

                        case "CF": // Correction formula
                            uniprotCF = ChemicalFormula.ParseFormula(line.Substring(5).Replace(" ", string.Empty));
                            break;

                        case "MM": // Monoisotopic mass difference. Might not precisely correspond to formula!
                            double thisMM;
                            if (!double.TryParse(line.Substring(5), out thisMM))
                                throw new PtmListLoaderException(line.Substring(5) + " is not a valid monoisotopic mass");
                            uniprotMM = thisMM;
                            break;

                        case "DR": // External database links!
                            var splitString = line.Substring(5).TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);
                            IList<string> val;
                            if (uniprotDR.TryGetValue(splitString[0], out val))
                                val.Add(splitString[1]);
                            else
                                uniprotDR.Add(splitString[0], new List<string> { splitString[1] });
                            break;

                        // NOW CUSTOM FIELDS:

                        case "NL": // Netural Losses. If field doesn't exist, single equal to 0
                            neutralLosses = new HashSet<double>(line.Substring(5).Split(new string[] { " or " }, StringSplitOptions.None).Select(b => double.Parse(b)));
                            break;

                        case "OM": // What masses are seen in histogram. If field doesn't exist, single equal to MM
                            massesObserved = new HashSet<double>(line.Substring(5).Split(new string[] { " or " }, StringSplitOptions.None).Select(b => double.Parse(b)));
                            break;

                        case "DI": // Masses of diagnostic ions
                            var nice = line.Substring(5).Split(new string[] { " or " }, StringSplitOptions.None);
                            if (!string.IsNullOrEmpty(nice[0]))
                                diagnosticIons = new HashSet<double>(nice.Select(b => double.Parse(b)));
                            else
                                diagnosticIons = null;
                            break;

                        case "MT": // Modification Type. If the field doesn't exist, set to the database name
                            modificationType = line.Substring(5);
                            break;

                        case "//":
                            // Not CROSSLNK, LIPID and MOD_RES is fine.
                            if ((uniprotFT == null || !uniprotFT.Equals("CROSSLNK")) && uniprotPP != null && uniprotTG != null && uniprotID != null)
                            {
                                if (ModificationWithLocation.modificationTypeCodes.TryGetValue(uniprotPP, out ModificationSites modSites))
                                {
                                    foreach (var singleTarget in uniprotTG)
                                    {
                                        string theMotif;
                                        if (aminoAcidCodes.TryGetValue(singleTarget, out char possibleMotifChar))
                                            theMotif = possibleMotifChar.ToString();
                                        else
                                            theMotif = singleTarget;
                                        if (ModificationMotif.TryGetMotif(theMotif, out ModificationMotif motif))
                                        {
                                            // Add the modification!
                                            if (!uniprotMM.HasValue)
                                            {
                                                // Return modification
                                                result = new ModificationWithLocation(uniprotID, uniprotAC, motif, modSites, uniprotDR, modificationType);
                                            }
                                            else
                                            {
                                                if (neutralLosses == null)
                                                    neutralLosses = new HashSet<double> { 0 };
                                                foreach (var neutralLoss in neutralLosses)
                                                {
                                                    if (uniprotCF == null)
                                                    {
                                                        // Return modification with mass
                                                        result = new ModificationWithMass(uniprotID, uniprotAC, motif, modSites, uniprotMM.Value, uniprotDR,
                                                            neutralLoss,
                                                            massesObserved ?? new HashSet<double> { uniprotMM.Value },
                                                            diagnosticIons,
                                                            modificationType);
                                                    }
                                                    else
                                                    {
                                                        // Return modification with complete information!
                                                        result = new ModificationWithMassAndCf(uniprotID, uniprotAC, motif, modSites, uniprotCF, uniprotMM.Value, uniprotDR,
                                                            neutralLoss,
                                                            massesObserved ?? new HashSet<double> { uniprotMM.Value },
                                                            diagnosticIons,
                                                            modificationType);
                                                    }
                                                }
                                            }
                                        }
                                        else
                                            throw new PtmListLoaderException("Could not get motif from " + singleTarget);
                                    }
                                }
                                else
                                    throw new PtmListLoaderException("Could not get modification site from " + uniprotPP);
                            }
                            break;
                    }
                }
            }
            return result;
        }

        #endregion Private Methods

    }
}