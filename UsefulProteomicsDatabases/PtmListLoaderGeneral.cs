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
    public static class PtmListLoaderGeneral
    {
        private static readonly Dictionary<string, char> aminoAcidCodes;

        static PtmListLoaderGeneral()
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

        public static IEnumerable<ModificationGeneral> ReadModsFromFile(string ptmListLocation)
        {
            return ReadModsFromFile(ptmListLocation, new Dictionary<string, int>()).OrderBy(b => b.Id);
        }

        /// <summary>
        /// Reads a list of modifications from a text file.
        /// </summary>
        /// <param name="ptmListLocation"></param>
        /// <returns></returns>
        public static IEnumerable<ModificationGeneral> ReadModsFromFile(string ptmListLocation, Dictionary<string, int> formalChargesDictionary)
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
                            yield return mod;
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
        public static IEnumerable<ModificationGeneral> ReadModsFromString(string storedModifications)
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
                            yield return mod;
                        modification_specification = new List<string>();
                    }
                }
            }
        }

        private static IEnumerable<ModificationGeneral> ReadMod(String ptmListLocation, List<string> specification, Dictionary<string, int> formalChargesDictionary)
        {
            string _Id = null;
            string _Accession = null;
            string _ModificationType = null;
            string _FeatureType = null;
            List<ModificationMotif> _Target = null;
            string _Position = null; //constructor will convert this to enum type
            ChemicalFormula _ChemicalFormula = null;
            double? _MonoisotopicMass = null;
            Dictionary<string, IList<string>> _DatabaseReference = new Dictionary<string, IList<string>>();
            Dictionary<string, IList<string>> _TaxonomicRange = new Dictionary<string, IList<string>>();
            List<string> _Keywords = null;
            Dictionary<DissociationType, List<double>> _NeutralLosses = new Dictionary<DissociationType, List<double>>();
            Dictionary<DissociationType, List<double>> _DiagnosticIons = new Dictionary<DissociationType, List<double>>();
            string _FileOrigin = ptmListLocation;

            foreach (string line in specification)
            {
                if (line.Length >= 2)
                {
                    string modKey = line.Substring(0, 2);
                    string modValue = null;
                    if (line.Length > 5)
                        try
                        {
                            modValue = line.Split('#')[0].Trim().Substring(5); //removes commented text
                        }
                        catch
                        {
                            //do nothing leave as null
                        }

                    switch (modKey)
                    {
                        case "ID": // Mandatory
                            _Id = modValue;
                            break;

                        case "AC": // Do not use! Only present in UniProt ptmlist
                            _Accession = modValue;
                            break;

                        case "FT": // Optional
                            _FeatureType = modValue;
                            break;

                        case "TG": // Which amino acid(s) or motifs is the modification on
                            string[] possibleMotifs = modValue.TrimEnd('.').Split(new string[] { " or " }, StringSplitOptions.None);
                            List<ModificationMotif> acceptableMotifs = new List<ModificationMotif>();
                            foreach (var singleTarget in possibleMotifs)
                            {
                                string theMotif;
                                if (aminoAcidCodes.TryGetValue(singleTarget, out char possibleMotifChar))
                                    theMotif = possibleMotifChar.ToString();
                                else
                                    theMotif = singleTarget;
                                if (ModificationMotif.TryGetMotif(theMotif, out ModificationMotif motif))
                                    acceptableMotifs.Add(motif);
                            }
                            _Target = acceptableMotifs.ToList();
                            break;

                        case "PP": // Terminus localization
                            _Position = modValue;
                            break;

                        case "CF": // Correction formula
                            _ChemicalFormula = ChemicalFormula.ParseFormula(modValue.Replace(" ", string.Empty));
                            break;

                        case "MM": // Monoisotopic mass difference. Might not precisely correspond to formula!
                            {
                                if (!double.TryParse(modValue, NumberStyles.Any, CultureInfo.InvariantCulture, out double thisMM))
                                    throw new MzLibException(line.Substring(5) + " is not a valid monoisotopic mass");
                                _MonoisotopicMass = thisMM;
                            }
                            break;

                        case "DR": // External database links!
                            {
                                var splitString = modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);
                                try
                                {
                                    _DatabaseReference.TryGetValue(splitString[0], out IList<string> val);
                                    val.Add(splitString[1]);
                                }
                                catch
                                {
                                    _DatabaseReference.Add(splitString[0], new List<string> { splitString[1] });
                                }
                            }
                            break;

                        case "TR": // External database links!
                            {
                                var splitString = modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);
                                try
                                {
                                    _TaxonomicRange.TryGetValue(splitString[0], out IList<string> val);
                                    val.Add(splitString[1]);
                                }
                                catch
                                {
                                    _TaxonomicRange.Add(splitString[0], new List<string> { splitString[1] });
                                }
                            }
                            break;

                        case "KW": // ; Separated keywords
                            {
                                _Keywords = new List<string>(modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None));
                            }
                            break;

                        // NOW CUSTOM FIELDS:

                        case "NL": // Netural Losses. If field doesn't exist, single equal to 0. these must all be on one line;
                            {
                                _NeutralLosses = DiagnosticIonsAndNeutralLosses(modValue);
                            }
                            break;

                        case "DI": // Masses of diagnostic ions. Might just be "DI"!!! If field doesn't exist, create an empty list!
                            {
                                _DiagnosticIons = DiagnosticIonsAndNeutralLosses(modValue);
                            }
                            break;

                        case "MT": // Modification Type. If the field doesn't exist, set to the database name
                            _ModificationType = modValue;
                            break;

                        case "//":

                            if (_Target.Count != 0)
                                foreach (ModificationMotif motif in _Target)
                                {
                                    yield return new ModificationGeneral(_Id, _Accession, _ModificationType, _FeatureType, motif, _Position, _ChemicalFormula, _MonoisotopicMass, _DatabaseReference, _TaxonomicRange, _Keywords, _NeutralLosses, _DiagnosticIons, _FileOrigin);
                                }
                            else
                            {
                                yield return new ModificationGeneral(_Id, _Accession, _ModificationType, _FeatureType, null, _Position, _ChemicalFormula, _MonoisotopicMass, _DatabaseReference, _TaxonomicRange, _Keywords, _NeutralLosses, _DiagnosticIons, _FileOrigin);
                            }
                            break;

                        default:
                            break;
                    }
                }
            }
        }

        public static DissociationType? ModDissociationType(string modType)
        {
            switch (modType)
            {
                case "Any":
                    return DissociationType.AnyActivationType;

                case "CID":
                    return DissociationType.CID;

                case "MPD":
                    return DissociationType.MPD;

                case "ECD":
                    return DissociationType.ECD;

                case "PQD":
                    return DissociationType.PQD;

                case "ETD":
                    return DissociationType.ETD;

                case "HCD":
                    return DissociationType.HCD;

                case "SA":
                    return DissociationType.SA;

                case "PTR":
                    return DissociationType.PTR;

                case "NETD":
                    return DissociationType.NETD;

                case "NPTR":
                    return DissociationType.NPTR;

                case "ISCID":
                    return DissociationType.ISCID;

                default:
                    return null;
            }
        }

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
                throw new MzLibException("woah");
            }

            return dAndNDictionary;
        }
    }
}