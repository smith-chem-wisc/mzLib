#nullable enable
using Omics.Modifications;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using Chemistry;
using Transcriptomics;
using System.Reflection;
using System.Diagnostics;
using CsvHelper;
using System.Globalization;

namespace UsefulProteomicsDatabases.Transcriptomics;

public static class ModomicsLoader
{
    // Json gives us most of the information 
    private static string _modomicsJsonPath = "UsefulProteomicsDatabases.Resources.modomics.json";

    // Csv tells us if the chemical formula is for a nucleotide or nucleoside. 
    private static string _modomicsCsvPath = "UsefulProteomicsDatabases.Resources.modomics.csv";
    internal static List<Modification>? ModomicsModifications { get; private set; }

    public static List<Modification> LoadModomics()
    {
        // Cache results
        if (ModomicsModifications != null)
            return ModomicsModifications;
        ModomicsModifications = new();

        // Load embedded resource
        var info = Assembly.GetExecutingAssembly().GetName().Name;
        var modomicsJsonStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{info}.Resources.modomicsmods.json");

        if (modomicsJsonStream is null)
            throw new FileNotFoundException("Could not find embedded resource", _modomicsJsonPath);
        var jsonDict = JsonSerializer.Deserialize<Dictionary<int, JsonElement>>(modomicsJsonStream)!;

        var modomicsCsvStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{info}.Resources.modomicsmods.csv");
        if (modomicsCsvStream is null)
            throw new FileNotFoundException("Could not find embedded resource", _modomicsCsvPath);

        // Parse CSV to get MoietyType (nucleotide vs nucleoside) mapping
        var moietyTypeByShortName = new Dictionary<string, string>();
        using (var reader = new StreamReader(modomicsCsvStream))
        using (var csv = new CsvReader(reader, CultureInfo.InvariantCulture))
        {
            csv.Read();
            csv.ReadHeader();
            while (csv.Read())
            {
                var shortName = csv.GetField("Short Name");
                var moietyType = csv.GetField("Moiety type");
                if (!string.IsNullOrWhiteSpace(shortName) && !string.IsNullOrWhiteSpace(moietyType))
                    moietyTypeByShortName[shortName] = moietyType;
            }
        }

        // Parse each entry into a DTO from the json
        var dtos = new List<ModomicsDto>();
        foreach (var kvp in jsonDict)
        {
            var modEntry = kvp.Value;

            // Declare each property, find the value, and then build the dto
            string abbrev = modEntry.GetProperty("abbrev").GetString()!;
            string formula = modEntry.GetProperty("formula").GetString()!.Replace("+", string.Empty);
            string name = modEntry.GetProperty("name").GetString()!;
            string? lcElutionComment = modEntry.TryGetProperty("lc_elution_comment", out var elutionComment) ? elutionComment.GetString() : null;
            string? lcElutionTime = modEntry.TryGetProperty("lc_elution_time", out var elutionTime) ? elutionTime.GetString() : null;
            double massAvg = modEntry.TryGetProperty("mass_avg", out var massAvgElem) && massAvgElem.ValueKind == JsonValueKind.Number
                ? massAvgElem.TryGetDouble(out var parsedMass) ? parsedMass : 0
                : 0;
            double massMonoiso = modEntry.TryGetProperty("mass_monoiso", out var massMonoElem) && massMonoElem.ValueKind == JsonValueKind.Number
                ? massMonoElem.TryGetDouble(out var parsedMass2) ? parsedMass2 : 0
                : 0;
            double massProt = modEntry.TryGetProperty("mass_prot", out var massProtElem) && massProtElem.ValueKind == JsonValueKind.Number
                ? massProtElem.TryGetDouble(out var parsedMass3) ? parsedMass3 : 0
                : 0;
            string? productIons = modEntry.TryGetProperty("product_ions", out var productIonsElem) ? productIonsElem.GetString()! : null;
            List<string?> referenceMoiety = modEntry.GetProperty("reference_moiety").EnumerateArray().Select(x => x.GetString()).Where(x => x != null).ToList();
            string shortName = modEntry.GetProperty("short_name").GetString()!;
            string? smile = modEntry.TryGetProperty("smile", out var smileElem) ? smileElem.GetString() : null;

            var dto = new ModomicsDto
            {
                Id = kvp.Key,
                Abbrev = abbrev,
                Formula = formula,
                LcElutionComment = lcElutionComment,
                LcElutionTime = lcElutionTime,
                MassAvg = massAvg,
                MassMonoiso = massMonoiso,
                MassProt = massProt,
                Name = name,
                ProductIons = productIons,
                ReferenceMoiety = referenceMoiety,
                ShortName = shortName,
                Smile = smile,
                MoietyType = moietyTypeByShortName.GetValueOrDefault(shortName)
            };

            dtos.Add(dto);
        }

        ModomicsModifications = dtos
            .SelectMany(dto => dto.ToModification())
            .ToList();

        return ModomicsModifications;
    }

    internal static IEnumerable<Modification> ToModification(this ModomicsDto dto)
    {
        string localizationRestriction = "Anywhere.";
        string modificationType = "Modomics";

        // Stand in mods for unknown residues, ignore them
        if (dto.Name.Contains("unknown", System.StringComparison.InvariantCultureIgnoreCase))
            yield break;

        // Caps and other weird mods. TODO: Handle those that are appropriate to do so
        if (dto.ReferenceMoiety.Count != 1)
            Debugger.Break();

        foreach (var refMoiety in dto.ReferenceMoiety)
        {
            // Filters out X as a moiety which occurs when mods are aimed at purine or pyrimidines TODO: Handle these appropriately
            if (!Nucleotide.AllKnownResidues.TryGetValue(refMoiety, out var baseNuc))
                continue;

            // Subtract nucleoside formula/mass
            var fullFormula = ChemicalFormula.ParseFormula(dto.Formula);
            ChemicalFormula modFormula = new ChemicalFormula(fullFormula);

            if (dto.MoietyType == "nucleoside")
            {
                modFormula.Remove(baseNuc.NucleosideChemicalFormula);
            }
            else if (dto.MoietyType == "nucleotide")
            {
                var sugarToRemove = baseNuc.NucleosideChemicalFormula - baseNuc.BaseChemicalFormula - Constants.WaterChemicalFormula;
                modFormula.Remove(sugarToRemove);
            }
            else if (dto.MoietyType == "base")
            {
                // These are unique bases due to modifications, will need to create a new nucleotide or create an x->This mod for each nucleotide. 
                yield break; 
            }
            else
            {
                // fallback or throw
                modFormula.Remove(baseNuc.NucleosideChemicalFormula);
            }
            //if (dto.Name.Contains("cap")) // Caps are in the database as mod + sugar
            //{
            //    var sugarToRemove = baseNuc.NucleosideChemicalFormula - baseNuc.BaseChemicalFormula - Constants.WaterChemicalFormula;
            //    modFormula.Remove(sugarToRemove);
            //}
            //else // The rest are in the database as full nucleotides. 
            //{
            //    modFormula.Remove(baseNuc.NucleosideChemicalFormula);
            //}

            if (!ModificationMotif.TryGetMotif(refMoiety, out var motif))
                continue;

            yield return new Modification(
                _originalId: dto.ShortName,
                _modificationType: modificationType,
                _target: motif,
                _locationRestriction: localizationRestriction,
                _chemicalFormula: modFormula,
                _monoisotopicMass: modFormula.MonoisotopicMass,
                _databaseReference: new Dictionary<string, IList<string>> { { "Modomics", new List<string> { dto.ShortName, dto.Name } } },
                _keywords: new List<string> { dto.Abbrev, dto.ShortName, dto.Name }
            );
        }
    }
}