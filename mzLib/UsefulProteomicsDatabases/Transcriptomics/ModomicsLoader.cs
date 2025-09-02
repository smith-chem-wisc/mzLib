#nullable enable
using Omics.Modifications;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using Chemistry;
using Transcriptomics;
using System.Reflection;

namespace UsefulProteomicsDatabases.Transcriptomics;

public static class ModomicsLoader
{
    private static string _modomicsResourceName = "UsefulProteomicsDatabases.Transcriptomics.modomics.json";
    public static List<Modification>? ModomicsModifications { get; private set; }

    public static List<Modification> LoadModomics()
    {
        // Cache results
        if (ModomicsModifications != null)
            return ModomicsModifications;
        ModomicsModifications = new();

        // Load embedded resource
        var info = Assembly.GetExecutingAssembly().GetName().Name;
        var modomicsStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{info}.Transcriptomics.modomics.json");

        if (modomicsStream is null)
            throw new FileNotFoundException("Could not find embedded resource", _modomicsResourceName);
        var modsDict = JsonSerializer.Deserialize<Dictionary<int, JsonElement>>(modomicsStream)!;

        // Parse each entry into a DTO
        var dtos = new List<ModomicsDto>();
        foreach (var kvp in modsDict)
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
                Smile = smile
            };

            dtos.Add(dto);
        }

        var dtosToParse = dtos
            .Where(p => !p.Name.Contains("unknown", System.StringComparison.InvariantCultureIgnoreCase));

        ModomicsModifications = dtosToParse
            .Select(dto => dto.ToModification())
            .Where(p => p is not null)
            .Cast<Modification>()
            .ToList();
        return ModomicsModifications;
    }


    public static Modification? ToModification(this ModomicsDto dto)
    {
        // Stand in mods for unknown residues, ignore them
        if (dto.Name.Contains("unknown", System.StringComparison.InvariantCultureIgnoreCase))
            return null!;

        // Caps and other weird mods. TODO: handle those that are appropriate to do so
        if (dto.ReferenceMoiety.Count > 1)
            return null!;



        // Defensive: Only use first reference moiety
        var refMoiety = dto.ReferenceMoiety?.FirstOrDefault();
        Nucleotide baseNuc = null;
        if (!string.IsNullOrEmpty(refMoiety) && Nucleotide.AllKnownResidues.TryGetValue(refMoiety, out var nuc))
            baseNuc = nuc;

        // Subtract nucleoside formula/mass if possible
        ChemicalFormula modFormula = null;
        double? modMonoMass = null;
        if (!string.IsNullOrWhiteSpace(dto.Formula) && baseNuc != null)
        {
            var fullFormula = ChemicalFormula.ParseFormula(dto.Formula);
            modFormula = new ChemicalFormula(fullFormula);
            modFormula.Remove(baseNuc.NucleosideChemicalFormula);
        }
        if (dto.MassMonoiso > 0 && baseNuc != null)
        {
            var baseMonoMass = baseNuc.NucleosideChemicalFormula.MonoisotopicMass;
            modMonoMass = dto.MassMonoiso - baseMonoMass;
        }

        ModificationMotif motif = null;
        if (!string.IsNullOrEmpty(refMoiety) && ModificationMotif.TryGetMotif(refMoiety, out var m))
            motif = m;

        return new Modification(
            _originalId: dto.ShortName,
            _modificationType: "RNA",
            _featureType: "MODOMICS",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: modFormula,
            _monoisotopicMass: modMonoMass,
            _databaseReference: new Dictionary<string, IList<string>> { { "MODOMICS", new List<string> { dto.ShortName } } },
            _keywords: new List<string> { dto.Abbrev, dto.Name }
        );
    }
}