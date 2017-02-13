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
            aminoAcidCodes = new Dictionary<string, char>();
            aminoAcidCodes.Add("Alanine", 'A');
            aminoAcidCodes.Add("Arginine", 'R');
            aminoAcidCodes.Add("Asparagine", 'N');
            aminoAcidCodes.Add("Aspartate", 'D');
            aminoAcidCodes.Add("Aspartic Acid", 'D');
            aminoAcidCodes.Add("Cysteine", 'C');
            aminoAcidCodes.Add("Glutamate", 'E');
            aminoAcidCodes.Add("Glutamic Acid", 'E');
            aminoAcidCodes.Add("Glutamine", 'Q');
            aminoAcidCodes.Add("Glycine", 'G');
            aminoAcidCodes.Add("Histidine", 'H');
            aminoAcidCodes.Add("Isoleucine", 'I');
            aminoAcidCodes.Add("Leucine", 'L');
            aminoAcidCodes.Add("Lysine", 'K');
            aminoAcidCodes.Add("Methionine", 'M');
            aminoAcidCodes.Add("Phenylalanine", 'F');
            aminoAcidCodes.Add("Proline", 'P');
            aminoAcidCodes.Add("Serine", 'S');
            aminoAcidCodes.Add("Threonine", 'T');
            aminoAcidCodes.Add("Tryptophan", 'W');
            aminoAcidCodes.Add("Tyrosine", 'Y');
            aminoAcidCodes.Add("Valine", 'V');
        }

        #endregion Public Constructors

        #region Public Methods

        public static IEnumerable<ModificationWithLocation> ReadMods(string ptmListLocation)
        {
            using (StreamReader uniprot_mods = new StreamReader(ptmListLocation))
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

                while (uniprot_mods.Peek() != -1)
                {
                    string line = uniprot_mods.ReadLine();
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

                            case "//":
                                // Only mod_res, not intrachain.
                                if ((uniprotFT == null || !uniprotFT.Equals("CROSSLNK")) && uniprotPP != null && uniprotTG != null && uniprotID != null)
                                {
                                    ModificationSites modSites;
                                    if (ModificationWithLocation.modificationTypeCodes.TryGetValue(uniprotPP, out modSites))
                                    {
                                        foreach (var singleTarget in uniprotTG)
                                        {
                                            string theMotif;
                                            char possibleMotifChar;
                                            if (aminoAcidCodes.TryGetValue(singleTarget, out possibleMotifChar))
                                                theMotif = possibleMotifChar.ToString();
                                            else
                                                theMotif = singleTarget;
                                            ModificationMotif motif;
                                            if (ModificationMotif.TryGetMotif(theMotif, out motif))
                                            {
                                                // Add the modification!
                                                if (!uniprotMM.HasValue)
                                                {
                                                    // Return modification
                                                    yield return new ModificationWithLocation(uniprotID, uniprotAC, motif, modSites, uniprotDR, Path.GetFileNameWithoutExtension(ptmListLocation));
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
                                                            yield return new ModificationWithMass(uniprotID, uniprotAC, motif, modSites, uniprotMM.Value, uniprotDR,
                                                                neutralLoss,
                                                                massesObserved == null ? new HashSet<double> { uniprotMM.Value } : massesObserved,
                                                                diagnosticIons,
                                                                Path.GetFileNameWithoutExtension(ptmListLocation));
                                                        }
                                                        else
                                                        {
                                                            // Return modification with complete information!
                                                            yield return new ModificationWithMassAndCf(uniprotID, uniprotAC, motif, modSites, uniprotCF, uniprotMM.Value, uniprotDR,
                                                                neutralLoss,
                                                                massesObserved == null ? new HashSet<double> { uniprotMM.Value } : massesObserved,
                                                                diagnosticIons,
                                                                Path.GetFileNameWithoutExtension(ptmListLocation));
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

                                uniprotID = null;
                                uniprotAC = null;
                                uniprotFT = null;
                                uniprotTG = null;
                                uniprotPP = null;
                                uniprotCF = null;
                                uniprotMM = null;
                                uniprotDR = new Dictionary<string, IList<string>>();

                                // Custom fields
                                neutralLosses = null;
                                massesObserved = null;
                                diagnosticIons = null;

                                break;
                        }
                    }
                }
            }
        }

        #endregion Public Methods

    }
}