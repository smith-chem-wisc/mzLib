using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics.Fragmentation;
using System.Globalization;
using System.Text.RegularExpressions;
using System.Xml.Serialization;

namespace Omics.Modifications.IO;
public static class ModificationLoader
{

    #region PTM List Loading (from ptmlist.txt format)

    private static readonly Dictionary<string, char> AminoAcidCodes = new()
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

    /// <summary>
    /// Reads a list of modifications from a text file.
    /// </summary>
    /// <param name="ptmListLocation"></param>
    /// <returns></returns>
    public static IEnumerable<Modification> ReadModsFromFile(string ptmListLocation, Dictionary<string, int> formalChargesDictionary, out List<(Modification, string)> filteredModificationsWithWarnings)
    {
        using (StreamReader uniprot_mods = new StreamReader(ptmListLocation))
        {
            return ReadModsFromFile(uniprot_mods, formalChargesDictionary, out filteredModificationsWithWarnings, ptmListLocation);
        }
    }
    
    /// <summary>
    /// Read a list of modifications from a text file.
    /// </summary>
    /// <param name="ptmListLocation"></param>
    /// <returns></returns>
    public static IEnumerable<Modification> ReadModsFromFile(string ptmListLocation, out List<(Modification, string)> filteredModificationsWithWarnings)
    {
        return ReadModsFromFile(ptmListLocation, new Dictionary<string, int>(), out filteredModificationsWithWarnings).OrderBy(b => b.IdWithMotif);
    }
    /// <summary>
    /// Reads a list of modifications from a stream reader.
    /// </summary>
    /// <param name="ptmListLocation"></param>
    /// <returns></returns>
    public static IEnumerable<Modification> ReadModsFromFile(StreamReader uniprot_mods, Dictionary<string, int> formalChargesDictionary, out List<(Modification, string)> filteredModificationsWithWarnings, string? fileLocation = null)
    {
        List<Modification> acceptedModifications = new List<Modification>();
        filteredModificationsWithWarnings = new List<(Modification filteredMod, string warningString)>();
        List<string> modification_specification = new List<string>();

        //This block will read one complete modification entry at a time until the EOF is reached.
        while (uniprot_mods.Peek() != -1) //The Peek method returns an integer value in order to determine whether the end of the file, or another error has occurred.
        {
            string line = uniprot_mods.ReadLine();
            modification_specification.Add(line);
            if (line.StartsWith("//"))
            {
                foreach (var mod in ReadMod(fileLocation, modification_specification, formalChargesDictionary))
                {
                    // Filter out the modifications that don't meet validation
                    if (mod.ValidModification)
                    {
                        acceptedModifications.Add(mod);
                    }
                    else
                    {
                        filteredModificationsWithWarnings.Add((mod, mod.ModificationErrorsToString()));
                    }
                }
                modification_specification = new List<string>();
            }
        }

        return acceptedModifications;
    }

    /// <summary>
    /// Reads a list of modifications from a string representation of a ptmlist text file.
    /// </summary>
    /// <param name="storedModifications"></param>
    /// <returns></returns>
    public static IEnumerable<Modification> ReadModsFromString(string storedModifications, out List<(Modification, string)> filteredModificationsWithWarnings)
    {
        List<Modification> acceptedModifications = new List<Modification>();
        filteredModificationsWithWarnings = new List<(Modification filteredMod, string warningString)>();
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
                        // Filter out the modifications that don't meet validation
                        if (mod.ValidModification)
                        {
                            acceptedModifications.Add(mod);
                        }
                        else
                        {
                            filteredModificationsWithWarnings.Add((mod, mod.ModificationErrorsToString()));
                        }
                    }
                    modification_specification = new List<string>();
                }
            }
        }
        return acceptedModifications;
    }

    /// <summary>
    /// Parse modification from string representation
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
        HashSet<ProductType>? _backboneProductTypes = null;
        BaseLossBehavior _baseLossType = BaseLossBehavior.Default;
        ChemicalFormula _baseLossModificationFormula = null;

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

                switch (modKey)
                {
                    case "ID": // Mandatory
                        _id = modValue;
                        break;

                    case "AC": // Do not use! Only present in UniProt ptmlist
                        _accession = modValue;
                        _modificationType = "UniProt";
                        break;

                    case "FT": // Optional
                        _featureType = modValue;
                        break;

                    case "TG": // Which amino acid(s) or motifs is the modification on
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
                        break;

                    case "PP": // Terminus localization
                        _locationRestriction = modValue;
                        break;

                    case "CF": // Correction formula
                        _chemicalFormula = ChemicalFormula.ParseFormula(modValue.Replace(" ", string.Empty));
                        break;

                    case "MM": // Monoisotopic mass difference. Might not precisely correspond to formula!
                        if (!double.TryParse(modValue, NumberStyles.Any, CultureInfo.InvariantCulture, out double thisMM))
                        {
                            throw new MzLibException(line.Substring(5) + " is not a valid monoisotopic mass");
                        }
                        _monoisotopicMass = thisMM;
                        break;

                    case "DR": // External database links!
                        var splitStringDR = modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);
                        if (_databaseReference == null)
                        {
                            _databaseReference = new Dictionary<string, IList<string>>();
                            _databaseReference.Add(splitStringDR[0], new List<string> { splitStringDR[1] });
                        }
                        else if (_databaseReference.TryGetValue(splitStringDR[0], out IList<string> val))
                        {
                            val.Add(splitStringDR[1]);
                        }
                        else
                        {
                            _databaseReference.Add(splitStringDR[0], new List<string> { splitStringDR[1] });
                        }
                        break;

                    case "TR": // External database links!
                        var splitStringTR = modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None);

                        if (_taxonomicRange == null)
                        {
                            _taxonomicRange = new Dictionary<string, IList<string>>();
                            _taxonomicRange.Add(splitStringTR[0], new List<string> { splitStringTR[1] });
                        }
                        else if (_taxonomicRange.TryGetValue(splitStringTR[0], out IList<string> val))
                        {
                            val.Add(splitStringTR[1]);
                        }
                        else
                        {
                            _taxonomicRange.Add(splitStringTR[0], new List<string> { splitStringTR[1] });
                        }
                        break;

                    case "KW": // ; Separated keywords
                        _keywords = new List<string>(modValue.TrimEnd('.').Split(new string[] { "; " }, StringSplitOptions.None));
                        break;

                    // NOW CUSTOM FIELDS:

                    case "NL": // Netural Losses. when field doesn't exist, single equal to 0. these must all be on one line;
                        if (_neutralLosses.IsNullOrEmpty())
                        {
                            _neutralLosses = new Dictionary<DissociationType, List<double>>();
                        }
                        _neutralLosses = DiagnosticIonsAndNeutralLosses(modValue, _neutralLosses);
                        break;

                    case "DI": // Masses of diagnostic ions. Might just be "DI"!!! If field doesn't exist, create an empty list!
                        if (_diagnosticIons.IsNullOrEmpty())
                        {
                            _diagnosticIons = new Dictionary<DissociationType, List<double>>();
                        }
                        _diagnosticIons = DiagnosticIonsAndNeutralLosses(modValue, _diagnosticIons);
                        break;

                    case "MT": // Modification Type. If the field doesn't exist, set to the database name
                        _modificationType = modValue;
                        break;

                    case "BM": // Backbone Modification (affects fragment masses)
                        if (string.IsNullOrWhiteSpace(modValue))
                            break;

                        // Parse: "b,c,d,x,y,z" - the fragments this mod can be added to 
                        string[] parts = modValue.Split(':');
                        string[] fragmentTypeStrings = parts[0].Split(',', StringSplitOptions.RemoveEmptyEntries);

                        // Parse fragment types
                        _backboneProductTypes = new HashSet<ProductType>();
                        foreach (string typeString in fragmentTypeStrings)
                        {
                            string trimmed = typeString.Trim();
                            if (Enum.TryParse(trimmed, true, out ProductType productType))
                            {
                                _backboneProductTypes.Add(productType);

                                foreach (var prodType in productType.GetFragmentFamilyMembers())
                                {
                                    _backboneProductTypes.Add(prodType);
                                }
                            }
                        }
                        break;

                    case "BL": // Base Loss behavior
                        if (string.IsNullOrWhiteSpace(modValue))
                            break;
                        // Parse optional formula
                        if (modValue.Contains(':'))
                        {
                            string formulaString = modValue.Substring(modValue.IndexOf(':') + 1).Trim();
                            try
                            {
                                _baseLossModificationFormula = ChemicalFormula.ParseFormula(formulaString.Replace(" ", string.Empty));
                            }
                            catch (Exception ex)
                            {
                                throw new MzLibException($"Invalid base loss formula '{formulaString}' for {_id}: {ex.Message}");
                            }
                        }

                        if (modValue.Equals("Suppressed", StringComparison.OrdinalIgnoreCase) ||
                            modValue.StartsWith("Suppressed:", StringComparison.OrdinalIgnoreCase))
                        {
                            _baseLossType = BaseLossBehavior.Suppressed;
                        }
                        else if (modValue.Equals("Modified", StringComparison.OrdinalIgnoreCase) ||
                                 modValue.StartsWith("Modified:", StringComparison.OrdinalIgnoreCase))
                        {
                            _baseLossType = BaseLossBehavior.Modified;
                            if (_baseLossModificationFormula == null)
                            {
                                // If the modification is a modified base loss but no formula is provided, default to the chemical formula of the modification
                                _baseLossModificationFormula = _chemicalFormula;
                            }
                        }
                        break;

                    case "//":
                        if (_target == null || _target.Count == 0) //This happens for FT=CROSSLINK modifications. We ignore these for now.
                        {
                            _target = new List<ModificationMotif> { null };
                        }
                        foreach (ModificationMotif motif in _target)
                        {
                            bool useChemFormulaForMM = _monoisotopicMass == null && _chemicalFormula != null;
                            bool adjustWithFormalCharge = _monoisotopicMass != null && _databaseReference != null;
                            if (useChemFormulaForMM)
                            {
                                _monoisotopicMass = _chemicalFormula.MonoisotopicMass;
                            }
                            if (adjustWithFormalCharge)
                            {
                                _monoisotopicMass = AdjustMonoIsotopicMassForFormalCharge(_monoisotopicMass, _chemicalFormula, _databaseReference, formalChargesDictionary);
                            }

                            if (_backboneProductTypes is { Count: > 0 })
                            {
                                yield return new BackboneModification(_id, _accession, _modificationType, _featureType, motif, _locationRestriction, _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin, _backboneProductTypes.ToArray());
                            }
                            else if (_baseLossType != BaseLossBehavior.Default || _baseLossModificationFormula != null)
                            {
                                yield return new BaseModification(_id, _accession, _modificationType, _featureType, motif, _locationRestriction, _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin, _baseLossType, _baseLossModificationFormula);
                            }
                            else
                            {

                                yield return new Modification(_id, _accession, _modificationType, _featureType, motif, _locationRestriction, _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin);
                            }
                        }
                        break;

                    default:
                        break;
                }
            }
        }
    }

    /// <summary>
    /// Subtract the mass of a proton for every formal charge on a modification.
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
    /// Parse diagnostic ion and neutral loss strings
    /// </summary>
    /// <param name="oneEntry"></param>
    /// <returns></returns>
    public static Dictionary<DissociationType, List<double>> DiagnosticIonsAndNeutralLosses(string oneEntry, Dictionary<DissociationType, List<double>> dAndNDictionary)
    {
        try
        {
            string[] nlOrDiEntries = oneEntry.Split(new string[] { " or " }, StringSplitOptions.None);
            foreach (string nlOrDiEntry in nlOrDiEntries)
            {
                string[] entryKeyValue = nlOrDiEntry.Split(new string[] { ":" }, StringSplitOptions.RemoveEmptyEntries);
                if (entryKeyValue.Length == 1) // assume there is no dissociation type listed and only formula or mass are supplied
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
                    if (val != null)
                    {
                        dAndNDictionary[DissociationType.AnyActivationType].Add(mm);
                    }
                    else
                    {
                        dAndNDictionary.Add(DissociationType.AnyActivationType, new List<double> { mm });
                    }
                }
                else if (entryKeyValue.Length == 2)  // an entry with two values is assumed to have a dissociation type and a neutral loss formula or mass
                {
                    DissociationType? dt = Enum.TryParse(entryKeyValue[0], true, out DissociationType parsedDt) ? parsedDt : null;

                    if (dt is null && entryKeyValue[0].Equals("MPD", StringComparison.InvariantCultureIgnoreCase))
                        dt = DissociationType.IRMPD;

                    if (dt != null)
                    {
                        //try // see if dictionary already contains key AnyActivationType
                        //{
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
                        if (val != null)
                        {
                            dAndNDictionary[(DissociationType)dt].Add(mm);
                        }
                        else
                        {
                            dAndNDictionary.Add((DissociationType)dt, new List<double> { mm });
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

    #endregion

    #region Unimod Loading (XML)

    private static readonly Dictionary<string, string> DictOfElements = new Dictionary<string, string>
        {
            {"2H", "H{2}" },
            {"13C", "C{13}" },
            {"18O", "O{18}" },
            {"15N", "N{15}" }
        };

    /// <summary>
    /// Read modifications from Unimod XML file
    /// </summary>
    public static IEnumerable<Modification> ReadModsFromUnimod(string unimodLocation)
    {
        using (FileStream stream = new FileStream(unimodLocation, FileMode.Open, FileAccess.Read, FileShare.Read))
        {
            return ReadModsFromUnimod(stream);
        }
    }

    /// <summary>
    /// Read modifications from Unimod XML stream
    /// </summary>
    public static IEnumerable<Modification> ReadModsFromUnimod(Stream unimodXml)
    {
        var unimodSerializer = new XmlSerializer(typeof(unimod_t));
        var deserialized = unimodSerializer.Deserialize(unimodXml) as unimod_t;
        return ReadModsFromUnimodDeserialized(deserialized);
    }

    private static IEnumerable<Modification> ReadModsFromUnimodDeserialized(
        unimod_t deserialized)
    {
        Dictionary<string, string> positionConversion = new Dictionary<string, string>()
            {
                { "Anywhere", "Anywhere."},
                { "AnyNterm", "Peptide N-terminal."},
                { "AnyCterm", "Peptide C-terminal."},
                { "ProteinNterm", "N-terminal."},
                { "ProteinCterm", "C-terminal."}
            };

        foreach (var mod in deserialized.modifications)
        {
            var id = mod.title;
            var ac = mod.record_id;
            ChemicalFormula cf = new ChemicalFormula();

            foreach (var el in mod.delta.element)
            {
                try
                {
                    cf.Add(el.symbol, int.Parse(el.number));
                }
                catch
                {
                    var tempCF = ChemicalFormula.ParseFormula(DictOfElements[el.symbol]);
                    tempCF.Multiply(int.Parse(el.number));
                    cf.Add(tempCF);
                }
            }

            foreach (var target in mod.specificity)
            {
                var tg = target.site;
                if (tg.Length > 1)
                    tg = "X";

                ModificationMotif.TryGetMotif(tg, out ModificationMotif motif);

                positionConversion.TryGetValue(target.position.ToString(), out string pos);

                Dictionary<string, IList<string>> dblinks = new Dictionary<string, IList<string>>
                    {
                        { "Unimod", new List<string>{ac.ToString() } },
                    };

                if (target.NeutralLoss == null)
                {
                    yield return new Modification(_originalId: id, _modificationType: "Unimod",
                        _target: motif, _locationRestriction: pos, _chemicalFormula: cf,
                        _databaseReference: dblinks);
                }
                else
                {
                    Dictionary<DissociationType, List<double>> neutralLosses = null;
                    foreach (var nl in target.NeutralLoss)
                    {
                        ChemicalFormula cfnl = new ChemicalFormula();
                        if (nl.mono_mass == 0)
                        {
                            if (neutralLosses == null)
                            {
                                neutralLosses = new Dictionary<DissociationType, List<double>>
                                    {
                                        { DissociationType.AnyActivationType, new List<double> { 0 }}
                                    };
                            }
                            else
                            {
                                if (neutralLosses.ContainsKey(DissociationType.AnyActivationType))
                                {
                                    if (!neutralLosses[DissociationType.AnyActivationType].Contains(0))
                                    {
                                        neutralLosses[DissociationType.AnyActivationType].Add(0);
                                    }
                                }
                            }
                        }
                        else
                        {
                            foreach (var el in nl.element)
                            {
                                try
                                {
                                    cfnl.Add(el.symbol, int.Parse(el.number));
                                }
                                catch
                                {
                                    var tempCF = ChemicalFormula.ParseFormula(DictOfElements[el.symbol]);
                                    tempCF.Multiply(int.Parse(el.number));
                                    cfnl.Add(tempCF);
                                }
                            }

                            if (neutralLosses == null)
                            {
                                neutralLosses = new Dictionary<DissociationType, List<double>>
                                    {
                                        { DissociationType.AnyActivationType, new List<double> { cfnl.MonoisotopicMass }}
                                    };
                            }
                            else
                            {
                                if (neutralLosses.ContainsKey(DissociationType.AnyActivationType))
                                {
                                    if (!neutralLosses[DissociationType.AnyActivationType].Contains(cfnl.MonoisotopicMass))
                                    {
                                        neutralLosses[DissociationType.AnyActivationType].Add(cfnl.MonoisotopicMass);
                                    }
                                }
                            }
                        }
                    }
                    yield return new Modification(_originalId: id, _target: motif,
                        _locationRestriction: "Anywhere.", _modificationType: "Unimod",
                        _chemicalFormula: cf, _databaseReference: dblinks, _neutralLosses: neutralLosses);
                }
            }
        }
    }

    #endregion

    #region PSI-MOD Loading

    /// <summary>
    /// Get formal charges dictionary from PSI-MOD XML
    /// </summary>
    public static Dictionary<string, int> GetFormalChargesDictionary(
        obo psiModDeserialized)
    {
        var modsWithFormalCharges = psiModDeserialized.Items
            .OfType<oboTerm>()
            .Where(b => b.xref_analog != null && b.xref_analog.Any(c => c.dbname.Equals("FormalCharge")));

        Regex digitsOnly = new(@"[^\d]");
        return modsWithFormalCharges.ToDictionary(
            b => "PSI-MOD; " + b.id,
            b => int.Parse(digitsOnly.Replace(
                b.xref_analog.First(c => c.dbname.Equals("FormalCharge")).name, "")));
    }

    /// <summary>
    /// Load PSI-MOD from XML file
    /// </summary>
    public static obo LoadPsiMod(string psimodLocation)
    {
        using FileStream stream = new FileStream(psimodLocation, FileMode.Open, FileAccess.Read, FileShare.Read);
        return LoadPsiMod(stream);
    }

    /// <summary>
    /// Load PSI-MOD from stream
    /// </summary>
    public static obo LoadPsiMod(Stream stream)
    {
        var psimodSerializer = new XmlSerializer(typeof(obo));
        return psimodSerializer.Deserialize(stream) as obo;
    }

    #endregion
}