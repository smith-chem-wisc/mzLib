using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    /// <summary>
    /// Class that contains methods for loading PTM lists from files or string values.
    /// </summary>
    public static class PtmListLoader
    {
        private static readonly Dictionary<string, char> AminoAcidCodes;

        static PtmListLoader()
        {
            AminoAcidCodes = new Dictionary<string, char>
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
                { "Pyrrolysine", 'O' },
                { "Selenocysteine", 'U' },
                { "Serine", 'S' },
                { "Threonine", 'T' },
                { "Tryptophan", 'W' },
                { "Tyrosine", 'Y' },
                { "Valine", 'V' }
            };
        }

        /// <summary>
        /// Reads a list of modifications from a text file.
        /// </summary>
        /// <param name="ptmListLocation"></param>
        /// <returns></returns>
        public static IEnumerable<Modification> ReadModsFromFile(string ptmListLocation)
        {
            return ReadModsFromFile(ptmListLocation, new Dictionary<string, int>()).OrderBy(b => b.IdWithMotif);
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
                        foreach (var mod in ReadMod(ptmListLocation, modification_specification, formalChargesDictionary))
                        {
                            // Filter out modifications that are invalid based on Proteomics.Modification.ValidModification
                            if (mod.ValidModification)
                            {
                                yield return mod;
                            }
                        }
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
                        foreach (var mod in ReadMod(null, modification_specification, new Dictionary<string, int>()))
                        {
                            yield return mod;
                        }
                        modification_specification = new List<string>();
                    }
                }
            }
        }

        /// <summary>
        /// Read information for a single modification based on string representation.
        /// </summary>
        /// <param name="ptmListLocation"></param>
        /// <param name="specification"></param>
        /// <param name="formalChargesDictionary"></param>
        /// <returns></returns>
        private static IEnumerable<Modification> ReadMod(string ptmListLocation, List<string> specification, Dictionary<string, int> formalChargesDictionary)
        {
            string _id = null;
            string _accession = null;
            string _modificationType = null;
            string _featureType = null;
            List<ModificationMotif> _target = null;
            string _locationRestriction = null; //constructor will convert this to enum type
            ChemicalFormula _chemicalFormula = null;
            double? _monoisotopicMass = null;
            Dictionary<string, IList<string>> _databaseReference = null;
            Dictionary<string, IList<string>> _taxonomicRange = null;
            List<string> _keywords = null;
            Dictionary<DissociationType, List<double>> _neutralLosses = null;
            Dictionary<DissociationType, List<double>> _diagnosticIons = null;
            string _fileOrigin = ptmListLocation;

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
                            modValue = line.Split('#')[0].Trim().Substring(5); //removes commented text
                        }
                        catch
                        {
                            //do nothing leave as null
                        }
                    }

                    if (modKey == "ID") // Mandatory
                    {
                        _id = modValue;
                    }

                    else if (modKey == "AC") // Do not use! Only present in UniProt ptmlist
                    {
                        _accession = modValue;
                        _modificationType = "UniProt";
                    }

                    else if (modKey == "FT") // Optional
                    {
                        _featureType = modValue;
                    }

                    else if (modKey == "TG") // Which amino acid(s) or motifs is the modification on
                    {
                        string[] possibleMotifs = modValue.TrimEnd('.').Split(new string[] { " or " }, StringSplitOptions.None);
                        List<ModificationMotif> acceptableMotifs = new List<ModificationMotif>();
                        foreach (var singleTarget in possibleMotifs)
                        {
                            string theMotif;
                            if (AminoAcidCodes.TryGetValue(singleTarget, out char possibleMotifChar))
                            {
                                theMotif = possibleMotifChar.ToString();
                            }
                            else
                            {
                                theMotif = singleTarget;
                            }
                            if (ModificationMotif.TryGetMotif(theMotif, out ModificationMotif motif))
                            {
                                acceptableMotifs.Add(motif);
                            }
                        }
                        _target = acceptableMotifs.ToList();
                    }

                    else if (modKey == "PP") // Terminus localization
                    {
                        _locationRestriction = modValue;
                    }

                    else if (modKey == "CF") // Correction formula
                    {
                        _chemicalFormula = ChemicalFormula.ParseFormula(modValue.Replace(" ", string.Empty));
                    }

                    else if (modKey == "MM") // Monoisotopic mass difference. Might not precisely correspond to formula!
                    {
                        if (!double.TryParse(modValue, NumberStyles.Any, CultureInfo.InvariantCulture, out double thisMM))
                        {
                            throw new MzLibException(line.Substring(5) + " is not a valid monoisotopic mass");
                        }
                        _monoisotopicMass = thisMM;
                    }

                    else if (modKey == "DR") // External database links!
                    {
                        var splitString = modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);
                        try
                        {
                            _databaseReference.TryGetValue(splitString[0], out IList<string> val);
                            val.Add(splitString[1]);
                        }
                        catch
                        {
                            if (_databaseReference == null)
                            {
                                _databaseReference = new Dictionary<string, IList<string>>();
                                _databaseReference.Add(splitString[0], new List<string> { splitString[1] });
                            }
                            else
                            {
                                _databaseReference.Add(splitString[0], new List<string> { splitString[1] });
                            }
                        }
                    }

                    else if (modKey == "TR") // External database links!
                    {
                        var splitString = modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);
                        try
                        {
                            _taxonomicRange.TryGetValue(splitString[0], out IList<string> val);
                            val.Add(splitString[1]);
                        }
                        catch
                        {
                            if (_taxonomicRange == null)
                            {
                                _taxonomicRange = new Dictionary<string, IList<string>>();
                                _taxonomicRange.Add(splitString[0], new List<string> { splitString[1] });
                            }
                            else
                            {
                                _taxonomicRange.Add(splitString[0], new List<string> { splitString[1] });
                            }
                        }
                    }

                    else if (modKey == "KW") // ; Separated keywords
                    {
                        _keywords = new List<string>(modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None));
                    }

                    // NOW CUSTOM FIELDS:
                    else if (modKey == "NL") // Netural Losses. when field doesn't exist, single equal to 0. these must all be on one line;
                    {
                        _neutralLosses = DiagnosticIonsAndNeutralLosses(modValue);
                    }

                    else if (modKey == "DI") // Masses of diagnostic ions. Might just be "DI"!!! If field doesn't exist, create an empty list!
                    {
                        _diagnosticIons = DiagnosticIonsAndNeutralLosses(modValue);
                    }

                    else if (modKey == "MT") // Modification Type. If the field doesn't exist, set to the database name
                    {
                        _modificationType = modValue;
                    }

                    else if (modKey == "//")
                    {
                        if (_target == null || _target.Count == 0) //This happens for FT=CROSSLINK modifications. We ignore these for now.
                        {
                            _target = new List<ModificationMotif> { null };
                        }
                        foreach (ModificationMotif motif in _target)
                        {
                            bool useChemicalFormulaForMM = _monoisotopicMass == null && _chemicalFormula != null;
                            bool adjustForFormalCharge = _monoisotopicMass != null && _databaseReference != null;

                            if (useChemicalFormulaForMM)
                            {
                                _monoisotopicMass = _chemicalFormula.MonoisotopicMass;
                            }
                            if (adjustForFormalCharge)
                            {
                                _monoisotopicMass = AdjustMonoIsotopicMassForFormalCharge(_monoisotopicMass, _chemicalFormula, _databaseReference, formalChargesDictionary);
                            }

                            yield return new Modification(_id, _accession, _modificationType, _featureType, motif, _locationRestriction, _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin);
                        }
                    }

                    else // tag not recognized
                    {
                        continue;
                    }
                }
            }
        }

        /// <summary>
        /// Adjust the monoisotopic mass by subtracting the mass of a proton for every formal charge.
        /// </summary>
        /// <param name="_monoisotopicMass"></param>
        /// <param name="_chemicalFormula"></param>
        /// <param name="_databaseReference"></param>
        /// <param name="formalChargesDictionary"></param>
        /// <returns></returns>
        private static double AdjustMonoIsotopicMassForFormalCharge(double? _monoisotopicMass, ChemicalFormula _chemicalFormula, Dictionary<string, IList<string>> _databaseReference, Dictionary<string, int> formalChargesDictionary)
        {
            foreach (var dbAndAccession in _databaseReference.SelectMany(b => b.Value.Select(c => b.Key + "; " + c)))
            {
                if (formalChargesDictionary.ContainsKey(dbAndAccession))
                {
                    if (_monoisotopicMass.HasValue)
                    {
                        _monoisotopicMass -= formalChargesDictionary[dbAndAccession] * Constants.ProtonMass;
                    }
                    if (_chemicalFormula != null)
                    {
                        _chemicalFormula.Remove(PeriodicTable.GetElement("H"), formalChargesDictionary[dbAndAccession]);
                    }
                    break;
                }
            }
            return (double)_monoisotopicMass;
        }

        /// <summary>
        /// Parse modification dissociation type string to enum.
        /// </summary>
        /// <param name="modType"></param>
        /// <returns></returns>
        public static DissociationType? ModDissociationType(string modType)
        {
            switch (modType)
            {
                case "Any":
                    return DissociationType.AnyActivationType;

                case "CID":
                    return DissociationType.CID;

                case "MPD":
                    return DissociationType.IRMPD;

                case "ECD":
                    return DissociationType.ECD;

                case "PQD":
                    return DissociationType.PQD;

                case "ETD":
                    return DissociationType.ETD;

                case "HCD":
                    return DissociationType.HCD;

                case "EThCD":
                    return DissociationType.EThcD;

                case "Custom":
                    return DissociationType.Custom;

                default:
                    return null;
            }
        }

        /// <summary>
        /// Parse diagnostic ion or neutral loss string representation
        /// </summary>
        /// <param name="oneEntry"></param>
        /// <returns></returns>
        public static Dictionary<DissociationType, List<double>> DiagnosticIonsAndNeutralLosses(string oneEntry)
        {
            Dictionary<DissociationType, List<double>> dAndNDictionary = new Dictionary<DissociationType, List<double>>();

            try
            {
                string[] nlOrDiEntries = oneEntry.Split(new string[] { " or " }, StringSplitOptions.None);
                foreach (string nlOrDiEntry in nlOrDiEntries)
                {
                    string[] entryKeyValue = nlOrDiEntry.Split(new string[] { ":" }, StringSplitOptions.RemoveEmptyEntries);
                    if (entryKeyValue.Length == 1) // assume there is no dissociation type listed and only formula or mass are supplied
                    {
                        try // see if dictionary already contains key AnyActivationType
                        {
                            double mm;
                            try
                            {
                                mm = ChemicalFormula.ParseFormula(entryKeyValue[0]).MonoisotopicMass; // turn chemical formula into monoisotopic mass
                            }
                            catch
                            {
                                mm = double.Parse(entryKeyValue[0], CultureInfo.InvariantCulture);
                            }

                            dAndNDictionary.TryGetValue(DissociationType.AnyActivationType, out List<double> val); // check the dictionary to see if AnyActivationType is already listed in the keys,
                            val.Add(mm); // no check for redundancy. could check for that in the future.
                        }
                        catch // dictionary doesn't contain key AnyActiviationType so we're gonna add it.
                        {
                            double mm;
                            try
                            {
                                mm = ChemicalFormula.ParseFormula(entryKeyValue[0]).MonoisotopicMass; // turn chemical formula into monoisotopic mass
                            }
                            catch
                            {
                                mm = double.Parse(entryKeyValue[0], CultureInfo.InvariantCulture);
                            }
                            dAndNDictionary.Add(DissociationType.AnyActivationType, new List<double>() { mm });
                        }
                    }
                    else if (entryKeyValue.Length == 2)  // an entry with two values is assumed to have a dissociation type and a neutral loss formula or mass
                    {
                        DissociationType? dt = ModDissociationType(entryKeyValue[0]);
                        if (dt != null)
                        {
                            try // see if dictionary already contains key AnyActivationType
                            {
                                double mm;
                                try
                                {
                                    mm = ChemicalFormula.ParseFormula(entryKeyValue[1]).MonoisotopicMass; // turn chemical formula into monoisotopic mass
                                }
                                catch
                                {
                                    mm = double.Parse(entryKeyValue[1], CultureInfo.InvariantCulture);
                                }

                                dAndNDictionary.TryGetValue((DissociationType)dt, out List<double> val); // check the dictionary to see if AnyActivationType is already listed in the keys,
                                val.Add(mm); // no check for redundancy. could check for that in the future.
                            }
                            catch // dictionary doesn't contain key AnyActiviationType so we're gonna add it.
                            {
                                double mm;
                                try
                                {
                                    mm = ChemicalFormula.ParseFormula(entryKeyValue[1]).MonoisotopicMass; // turn chemical formula into monoisotopic mass
                                }
                                catch
                                {
                                    mm = double.Parse(entryKeyValue[1], CultureInfo.InvariantCulture);
                                }
                                dAndNDictionary.Add((DissociationType)dt, new List<double>() { mm });
                            }
                        }
                        else
                        {
                            throw new MzLibException("neutral loss or diagnostic ion entry dissociation type is not parsable");
                        }
                    }
                    else
                    {
                        throw new MzLibException("your neutral loss or diagnostic ion is junk");
                    }
                }
            }
            catch
            {
                dAndNDictionary = null; // must have run into some junk
            }

            return dAndNDictionary;
        }
    }
}