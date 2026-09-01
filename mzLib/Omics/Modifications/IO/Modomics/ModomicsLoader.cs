using System.Globalization;
using System.Text.Json;
using Microsoft.VisualBasic.FileIO;

using Chemistry;
using MassSpectrometry;

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
    };

    // The ribose carried inside a nucleoside (nucleoside - bonded base), and the two hydrogens added when
    // a nucleoside fragments to its protonated free-base product ion (de-glycosylation + protonation).
    private static readonly ChemicalFormula NucleosideRiboseChemicalFormula = Formulas.AdenosineChemicalFormula - Formulas.AdenineBaseChemicalFormula;
    private static readonly ChemicalFormula TwoHydrogenChemicalFormula = ChemicalFormula.ParseFormula("H2");

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
            Formula = modEntry.GetProperty("formula").GetString() ?? string.Empty,
            LcElutionComment = GetStringOrNull(modEntry, "lc_elution_comment"),
            LcElutionTime = GetStringOrNull(modEntry, "lc_elution_time"),
            MassAvg = GetDoubleOrNull(modEntry, "mass_avg"),
            MassMonoiso = GetDoubleOrNull(modEntry, "mass_monoiso"),
            MassProt = GetDoubleOrNull(modEntry, "mass_prot"),
            Name = modEntry.GetProperty("name").GetString() ?? string.Empty,
            ProductIons = GetStringOrNull(modEntry, "product_ions"),
            ReferenceMoieties = modEntry.GetProperty("reference_moiety")
                .EnumerateArray()
                .Select(p => p.GetString())
                .OfType<string>()
                .ToList(),
            ShortName = shortName,
            Smile = GetStringOrNull(modEntry, "smile"),
            MoietyType = moietyTypeByShortName.GetValueOrDefault(shortName, string.Empty),
        };
    }

    private static string? GetStringOrNull(JsonElement element, string propertyName)
    {
        return element.TryGetProperty(propertyName, out var value) && value.ValueKind == JsonValueKind.String
            ? value.GetString()
            : null;
    }

    private static double? GetDoubleOrNull(JsonElement element, string propertyName)
    {
        return element.TryGetProperty(propertyName, out var value)
            && value.ValueKind == JsonValueKind.Number
            && value.TryGetDouble(out var mass)
            ? mass
            : null;
    }

    private static IEnumerable<ModomicsConversionOutcome> ConvertDto(ModomicsDto dto)
    {
        var neutralFormulaText = dto.Formula.Replace("+", string.Empty, StringComparison.Ordinal);

        if (string.IsNullOrWhiteSpace(neutralFormulaText))
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
            fullFormula = ChemicalFormula.ParseFormula(neutralFormulaText);
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

        if (!isTerminalCap && dto.ReferenceMoieties.Count > 1)
        {
            // Generic-"N" structures (NAD/CoA-linked 5' cofactors) reference every base and describe a
            // terminal cofactor rather than a residue mass shift.
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.TerminalStructureRequiresTerminusModel, dto.ShortName);
            yield break;
        }

        if (dto.ReferenceMoieties.Count == 0)
        {
            yield return NotRepresentable(dto, ModomicsRepresentationFailureReason.UnsupportedReferenceMoiety, "X");
            yield break;
        }

        var diagnosticIons = BuildDiagnosticIons(dto, fullFormula);

        foreach (var referenceMoiety in dto.ReferenceMoieties)
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
                Modification = CreateModification(dto, motif, isTerminalCap, modFormula, diagnosticIons),
            };
        }
    }

    private static Dictionary<DissociationType, List<double>>? BuildDiagnosticIons(ModomicsDto dto, ChemicalFormula fullFormula)
    {
        if (string.IsNullOrWhiteSpace(dto.ProductIons))
        {
            return null;
        }

        var ions = dto.ProductIons
            .Split('/', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .Select(token => double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out var mz) ? mz : (double?)null)
            .OfType<double>()
            .Distinct()
            .ToList();

        if (ions.Count == 0)
        {
            return null;
        }

        // MODOMICS publishes product ions as nominal integer m/z of protonated fragments from [M+H]+.
        // Proteomics diagnostic ions carry monoisotopic precision, so the primary base ion is upgraded to
        // its computed monoisotopic m/z whenever it validates against the published nominal value; ions we
        // cannot derive (secondary neutral-loss fragments, unmodified base ions, protonation-ambiguous
        // entries) keep their published nominal values.
        if (dto.MoietyType == "nucleoside")
        {
            var protonatedBaseIonMass = ComputeProtonatedBaseIonMass(fullFormula);
            if (protonatedBaseIonMass is not null)
            {
                var nominalBaseIon = (int)Math.Round(protonatedBaseIonMass.Value, MidpointRounding.AwayFromZero);
                var ionIndex = ions.FindIndex(mz => (int)Math.Round(mz, MidpointRounding.AwayFromZero) == nominalBaseIon);
                if (ionIndex >= 0)
                {
                    ions[ionIndex] = protonatedBaseIonMass.Value;
                }
            }
        }

        return new Dictionary<DissociationType, List<double>>
        {
            { DissociationType.AnyActivationType, ions },
        };
    }

    private static double? ComputeProtonatedBaseIonMass(ChemicalFormula fullFormula)
    {
        var baseIon = new ChemicalFormula(fullFormula);
        // The dominant nucleoside product ion is the protonated free base: remove the ribose carried in
        // the nucleoside, then add the two hydrogens (de-glycosylation + protonation).
        baseIon.Remove(NucleosideRiboseChemicalFormula);
        baseIon.Add(TwoHydrogenChemicalFormula);

        return string.IsNullOrEmpty(baseIon.Formula) ? null : baseIon.MonoisotopicMass;
    }

    private static Modification CreateModification(ModomicsDto dto, ModificationMotif motif, bool isTerminalCap, ChemicalFormula chemicalFormula, Dictionary<DissociationType, List<double>>? diagnosticIons)
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
            _keywords: new List<string> { dto.Abbrev, dto.ShortName, dto.Name },
            _diagnosticIons: diagnosticIons
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
                ReferenceMoieties = dto.ReferenceMoieties,
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
