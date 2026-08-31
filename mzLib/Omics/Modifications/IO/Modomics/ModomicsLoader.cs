using System.Globalization;
using System.Text.Json;
using Microsoft.VisualBasic.FileIO;

using Chemistry;

namespace Omics.Modifications.IO.Modomics;

/// <summary>
/// Loads MODOMICS RNA modifications from the embedded database and converts whole-nucleoside/nucleotide
/// formulas into mass shifts relative to the canonical base residues. Entries that cannot be represented
/// as a residue mass shift are reported rather than silently dropped.
/// </summary>
public static class ModomicsLoader
{
    private const string ModomicsJsonResourceName = "Resources.modomicsmods.json";
    private const string ModomicsCsvResourceName = "Resources.modomicsmods.csv";
    private static readonly object _cacheLock = new();
    private static ModomicsLoadResult? _cachedBaseResult;

    // Bonded base and neutral nucleoside formulas per reference moiety, from Chemistry.Formulas — the
    // shared ground truth for residue chemistry, also consumed by Transcriptomics.Nucleotide.
    private static readonly Dictionary<string, ReferenceMoietyDefinition> ReferenceMoieties = new(StringComparer.Ordinal)
    {
        { "A", new ReferenceMoietyDefinition(Formulas.AdenineBaseChemicalFormula, Formulas.AdenosineChemicalFormula) },
        { "C", new ReferenceMoietyDefinition(Formulas.CytosineBaseChemicalFormula, Formulas.CytidineChemicalFormula) },
        { "G", new ReferenceMoietyDefinition(Formulas.GuanineBaseChemicalFormula, Formulas.GuanosineChemicalFormula) },
        { "U", new ReferenceMoietyDefinition(Formulas.UracilBaseChemicalFormula, Formulas.UridineChemicalFormula) },
        { "Y", new ReferenceMoietyDefinition(Formulas.UracilBaseChemicalFormula, Formulas.UridineChemicalFormula) },
    };

    /// <summary>
    /// Loads MODOMICS RNA modifications.
    /// </summary>
    /// <param name="existingMetaMorpheusRnaMods">
    /// Optional curated MetaMorpheus RNA mods; when provided, chemistry-equivalent overlaps are recorded in
    /// <see cref="ModomicsLoadResult.DuplicateModifications"/>. Both entries are kept in the registry.
    /// </param>
    /// <returns>The structured load report.</returns>
    public static ModomicsLoadResult LoadModomics(IEnumerable<Modification>? existingMetaMorpheusRnaMods = null)
    {
        EnsureBaseCache();

        if (existingMetaMorpheusRnaMods is null)
        {
            return _cachedBaseResult!;
        }

        var duplicates = _cachedBaseResult!.LoadedModifications
            .Select(modomics => new
            {
                ModomicsModification = modomics,
                ExistingModifications = existingMetaMorpheusRnaMods
                    .Where(existing => AreEquivalent(existing, modomics))
                    .ToList(),
            })
            .Where(p => p.ExistingModifications.Count > 0)
            .Select(p => new ModomicsDuplicateEntry
            {
                ExistingModifications = p.ExistingModifications,
                ModomicsModification = p.ModomicsModification,
            })
            .ToList();

        return new ModomicsLoadResult
        {
            LoadedModifications = _cachedBaseResult.LoadedModifications,
            TerminalModifications = _cachedBaseResult.TerminalModifications,
            DuplicateModifications = duplicates,
            NotYetRepresentableEntries = _cachedBaseResult.NotYetRepresentableEntries,
        };
    }

    private static void EnsureBaseCache()
    {
        if (_cachedBaseResult is not null)
        {
            return;
        }

        lock (_cacheLock)
        {
            if (_cachedBaseResult is not null)
            {
                return;
            }

            _cachedBaseResult = BuildBaseResult();
        }
    }

    private static ModomicsLoadResult BuildBaseResult()
    {
        var assembly = typeof(ModomicsLoader).Assembly;
        var assemblyName = assembly.GetName().Name;
        using var modomicsJsonStream = assembly.GetManifestResourceStream($"{assemblyName}.{ModomicsJsonResourceName}")
            ?? throw new FileNotFoundException("Could not find embedded MODOMICS JSON resource", ModomicsJsonResourceName);
        using var modomicsCsvStream = assembly.GetManifestResourceStream($"{assemblyName}.{ModomicsCsvResourceName}")
            ?? throw new FileNotFoundException("Could not find embedded MODOMICS CSV resource", ModomicsCsvResourceName);

        var jsonDict = JsonSerializer.Deserialize<Dictionary<int, JsonElement>>(modomicsJsonStream)
            ?? throw new InvalidDataException("MODOMICS JSON resource could not be parsed.");

        var moietyTypeByShortName = ReadMoietyTypesByShortName(modomicsCsvStream);
        var loadedModifications = new List<Modification>();
        var terminalModifications = new List<Modification>();
        var notYetRepresentableEntries = new List<ModomicsNotYetRepresentableEntry>();
        var loadedById = new Dictionary<string, Modification>(StringComparer.Ordinal);

        foreach (var kvp in jsonDict)
        {
            var dto = BuildDto(kvp.Key, kvp.Value, moietyTypeByShortName);
            foreach (var outcome in ConvertDto(dto))
            {
                if (outcome.NotYetRepresentableEntry is not null)
                {
                    notYetRepresentableEntries.Add(outcome.NotYetRepresentableEntry);
                    continue;
                }

                if (outcome.Modification is null)
                {
                    continue;
                }

                if (loadedById.TryAdd(outcome.Modification.IdWithMotif, outcome.Modification))
                {
                    loadedModifications.Add(outcome.Modification);
                    if (outcome.Modification.LocationRestriction is "5'-terminal." or "Oligo 5'-terminal.")
                    {
                        terminalModifications.Add(outcome.Modification);
                    }
                }
            }
        }

        return new ModomicsLoadResult
        {
            LoadedModifications = loadedModifications,
            TerminalModifications = terminalModifications,
            NotYetRepresentableEntries = notYetRepresentableEntries,
        };
    }

    private static Dictionary<string, string> ReadMoietyTypesByShortName(Stream modomicsCsvStream)
    {
        var moietyTypeByShortName = new Dictionary<string, string>(StringComparer.Ordinal);
        using var parser = new TextFieldParser(modomicsCsvStream)
        {
            TextFieldType = FieldType.Delimited,
            HasFieldsEnclosedInQuotes = true,
        };
        parser.SetDelimiters(",");

        var header = parser.ReadFields() ?? [];
        var shortNameIndex = Array.FindIndex(header, p => string.Equals(p, "Short Name", StringComparison.Ordinal));
        var moietyTypeIndex = Array.FindIndex(header, p => string.Equals(p, "Moiety type", StringComparison.Ordinal));

        if (shortNameIndex < 0 || moietyTypeIndex < 0)
        {
            throw new InvalidDataException("MODOMICS CSV is missing required columns.");
        }

        while (!parser.EndOfData)
        {
            var fields = parser.ReadFields();
            if (fields is null || fields.Length <= Math.Max(shortNameIndex, moietyTypeIndex))
            {
                continue;
            }

            var shortName = fields[shortNameIndex];
            var moietyType = fields[moietyTypeIndex];
            if (!string.IsNullOrWhiteSpace(shortName) && !string.IsNullOrWhiteSpace(moietyType))
            {
                moietyTypeByShortName[shortName] = moietyType;
            }
        }

        return moietyTypeByShortName;
    }

    private static ModomicsDto BuildDto(int id, JsonElement modEntry, IReadOnlyDictionary<string, string> moietyTypeByShortName)
    {
        var shortName = modEntry.GetProperty("short_name").GetString() ?? string.Empty;
        return new ModomicsDto
        {
            Id = id,
            Abbrev = modEntry.GetProperty("abbrev").GetString() ?? string.Empty,
            Formula = (modEntry.GetProperty("formula").GetString() ?? string.Empty).Replace("+", string.Empty, StringComparison.Ordinal),
            Name = modEntry.GetProperty("name").GetString() ?? string.Empty,
            ReferenceMoiety = modEntry.GetProperty("reference_moiety")
                .EnumerateArray()
                .Select(p => p.GetString())
                .OfType<string>()
                .ToList(),
            ShortName = shortName,
            MoietyType = moietyTypeByShortName.GetValueOrDefault(shortName, string.Empty),
        };
    }

    private static IEnumerable<ModomicsConversionOutcome> ConvertDto(ModomicsDto dto)
    {
        if (string.IsNullOrWhiteSpace(dto.Formula))
        {
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.EmptyFormula);
            yield break;
        }

        if (dto.Name.Contains("unknown", StringComparison.InvariantCultureIgnoreCase))
        {
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.UnknownEntry);
            yield break;
        }

        ChemicalFormula? fullFormula = null;
        string? invalidFormulaMessage = null;
        try
        {
            fullFormula = ChemicalFormula.ParseFormula(dto.Formula);
        }
        catch (Exception ex)
        {
            invalidFormulaMessage = ex.Message;
        }

        if (fullFormula is null)
        {
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.InvalidFormula, invalidFormulaMessage);
            yield break;
        }

        var isTerminalCap = dto.Name.Contains(" cap", StringComparison.InvariantCultureIgnoreCase);

        if (dto.Name.Contains(" end", StringComparison.InvariantCultureIgnoreCase))
        {
            // Terminus-state descriptors ("5' hydroxyl end", "5' diphosphate end", ...) belong to the
            // five/three prime terminus model, not to a residue mass shift.
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.TerminalStructureRequiresTerminusModel, dto.ShortName);
            yield break;
        }

        if (!isTerminalCap && dto.ReferenceMoiety.Count > 1)
        {
            // Generic-"N" structures (NAD/CoA-linked 5' cofactors) reference every base and describe a
            // terminal cofactor rather than a residue mass shift.
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.TerminalStructureRequiresTerminusModel, dto.ShortName);
            yield break;
        }

        if (dto.ReferenceMoiety.Count == 0)
        {
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.UnsupportedReferenceMoiety, "X");
            yield break;
        }

        foreach (var referenceMoiety in dto.ReferenceMoiety)
        {
            if (!ReferenceMoieties.TryGetValue(referenceMoiety, out var moietyDefinition))
            {
                // "X" denotes generic purine/pyrimidine targets, which require per-nucleotide expansion
                // before they can be represented as a mass shift on a specific residue.
                yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.UnsupportedReferenceMoiety, referenceMoiety);
                continue;
            }

            if (!ModificationMotif.TryGetMotif(referenceMoiety, out var motif))
            {
                yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.InvalidTargetMotif, referenceMoiety);
                continue;
            }

            ChemicalFormula referenceFormulaToRemove;
            switch (dto.MoietyType)
            {
                case "nucleoside":
                    // MODOMICS nucleoside formulas describe the entire (neutral) nucleoside; the residue
                    // already carries the canonical nucleoside scaffold, so removing it yields the shift.
                    referenceFormulaToRemove = moietyDefinition.NucleosideChemicalFormula;
                    break;
                case "nucleotide" when isTerminalCap:
                    // Generic-"N" cap formulas embed a base-less ribose for the capped residue; removing only
                    // that ribose keeps the cap nucleoside and phosphate chain as the terminal mass shift.
                    referenceFormulaToRemove = moietyDefinition.NucleosideChemicalFormula - moietyDefinition.BaseChemicalFormula - Formulas.WaterChemicalFormula;
                    break;
                case "nucleotide":
                    // Real-base nucleotide formulas are published as inconsistent phosphate anions (e.g. pm1G
                    // is "C11H14N5O8P", two protons shy of the free acid), so an exact neutral shift cannot
                    // be derived reliably; the neutral nucleoside entries carry the same chemistry.
                    yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.NucleotideProtonationAmbiguous, dto.MoietyType);
                    continue;
                case "base":
                    // A modified base (e.g. preQ0) requires a new residue definition rather than a mass
                    // shift on an existing residue.
                    yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.BaseMoietyRequiresNewResidue, referenceMoiety);
                    continue;
                default:
                    yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.UnsupportedMoietyType, dto.MoietyType);
                    continue;
            }

            var modFormula = new ChemicalFormula(fullFormula);
            modFormula.Remove(referenceFormulaToRemove);

            if (string.IsNullOrEmpty(modFormula.Formula))
            {
                // Canonical nucleoside entries (e.g. "G on G") and isomerizations such as pseudouridine
                // impart no mass shift, so they cannot be represented as mass-difference modifications.
                yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.NoMassShift, referenceMoiety);
                continue;
            }

            yield return new ModomicsConversionOutcome
            {
                Modification = CreateModification(dto, motif, isTerminalCap, modFormula),
            };
        }
    }

    private static Modification CreateModification(ModomicsDto dto, ModificationMotif motif, bool isTerminalCap, ChemicalFormula chemicalFormula)
    {
        var idString = dto.Id.ToString(CultureInfo.InvariantCulture);

        return new Modification(
            _originalId: dto.Name,
            _accession: idString,
            _modificationType: isTerminalCap ? "5' Terminal Cap" : "Modomics",
            _target: motif,
            _locationRestriction: isTerminalCap ? "5'-terminal." : "Anywhere.",
            _chemicalFormula: chemicalFormula,
            _monoisotopicMass: chemicalFormula.MonoisotopicMass,
            _databaseReference: new Dictionary<string, IList<string>>
            {
                { "Modomics", new List<string> { dto.ShortName, dto.Name, idString } },
            },
            _keywords: new List<string> { dto.Abbrev, dto.ShortName, dto.Name }
        );
    }

    private static ModomicsConversionOutcome NotRepresentable(ModomicsDto dto, ModomicsRepresentationFailureReason reason, string? details = null)
    {
        return new ModomicsConversionOutcome
        {
            NotYetRepresentableEntry = new ModomicsNotYetRepresentableEntry
            {
                Id = dto.Id,
                ShortName = dto.ShortName,
                Name = dto.Name,
                Formula = dto.Formula,
                MoietyType = dto.MoietyType,
                ReferenceMoieties = dto.ReferenceMoiety,
                Reason = reason,
                Details = details,
            },
        };
    }

    private static bool AreEquivalent(Modification existingModification, Modification modomicsModification)
    {
        return existingModification.Target?.ToString() == modomicsModification.Target?.ToString()
               && existingModification.LocationRestriction == modomicsModification.LocationRestriction
               && existingModification.ChemicalFormula?.Equals(modomicsModification.ChemicalFormula) == true;
    }

    private sealed class ModomicsDto
    {
        public int Id { get; init; }
        public string Abbrev { get; init; } = string.Empty;
        public string Formula { get; init; } = string.Empty;
        public string Name { get; init; } = string.Empty;
        public List<string> ReferenceMoiety { get; init; } = [];
        public string ShortName { get; init; } = string.Empty;
        public string MoietyType { get; init; } = string.Empty;
    }

    private sealed class ModomicsConversionOutcome
    {
        public Modification? Modification { get; init; }
        public ModomicsNotYetRepresentableEntry? NotYetRepresentableEntry { get; init; }
    }

    private sealed class ReferenceMoietyDefinition(ChemicalFormula baseFormula, ChemicalFormula nucleosideFormula)
    {
        /// <summary>Bonded base formula (free base less one H), from Chemistry.Formulas.</summary>
        public ChemicalFormula BaseChemicalFormula { get; } = baseFormula;

        /// <summary>Neutral nucleoside formula (bonded base + ribose, condensed), from Chemistry.Formulas.</summary>
        public ChemicalFormula NucleosideChemicalFormula { get; } = nucleosideFormula;
    }
}
