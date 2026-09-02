using System.Globalization;
using System.Text.Json;
using System.Text.RegularExpressions;
using Microsoft.VisualBasic.FileIO;

using Chemistry;
using MassSpectrometry;
using MzLibUtil;

namespace Omics.Modifications.IO.Modomics;

/// <summary>
/// Loads MODOMICS RNA modifications from the embedded database and converts whole-nucleoside/nucleotide
/// formulas into mass shifts relative to the canonical base residues. Entries that cannot be represented
/// as a residue mass shift are reported rather than silently dropped.
/// Every published product ion is retained as an AnyActivationType diagnostic ion. For nucleoside
/// entries, the primary base cation (the largest published ion) additionally measures how much of the
/// modification sits on the base, which yields <see cref="BaseModification"/> semantics: nothing on the
/// base (a ribose methyl) gives Default; a partial or full portion gives Modified with that portion
/// departing with the base. Ions that measure no unique answer keep the plain representation.
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

    // The two hydrogens added when a nucleoside fragments to its protonated free-base product ion
    // (de-glycosylation + protonation): [base+H]+ = bonded base + 2H.
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

        var parsedProductIonCandidates = dto.ProductIons?
            .Split('/', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries)
            .Select(token => double.TryParse(token, NumberStyles.Float, CultureInfo.InvariantCulture, out var mz) ? mz : (double?)null)
            .OfType<double>()
            .Distinct()
            .ToList();

        // A present-but-empty or unparsable product_ions field yields no usable ions.
        var parsedProductIons = parsedProductIonCandidates is { Count: > 0 } ? parsedProductIonCandidates : null;

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

            var productIonInterpretation = InterpretProductIons(parsedProductIons, dto.MoietyType, moietyDefinition, modFormula);

            yield return new ModomicsConversionOutcome
            {
                Modification = CreateModification(dto, motif, isTerminalCap, modFormula, productIonInterpretation),
            };
        }
    }

    /// <summary>
    /// The loader's interpretation of a MODOMICS entry's published product ions: the ions to publish as
    /// diagnostic ions, and the base-loss semantics the primary ion measures.
    /// </summary>
    private sealed record ProductIonInterpretation(
        Dictionary<DissociationType, List<double>>? DiagnosticIons,
        BaseLossBehavior? BaseLossType,
        ChemicalFormula? BaseLossModification);

    /// <summary>
    /// Keeps every published ion as an AnyActivationType diagnostic ion. For nucleoside entries, the
    /// primary ion (the base cation) additionally measures how much of the modification sits on the
    /// base, which yields base-loss semantics and the primary ion's accurate monoisotopic m/z.
    /// </summary>
    private static ProductIonInterpretation? InterpretProductIons(
        List<double>? productIons, string moietyType, ReferenceMoietyDefinition moiety, ChemicalFormula modificationFormula)
    {
        if (productIons is null)
        {
            return null;
        }

        // Every published ion becomes a diagnostic ion.
        var ions = productIons.ToList();
        var diagnosticIons = new Dictionary<DissociationType, List<double>>
        {
            { DissociationType.AnyActivationType, ions },
        };

        var measured = moietyType == "nucleoside"
            ? DeriveBaseLossModification(productIons, (moiety.BaseChemicalFormula + TwoHydrogenChemicalFormula).MonoisotopicMass, modificationFormula)
            : null;

        if (measured is null)
        {
            // Caps, and nucleosides whose ions do not measure a unique base-localized portion: keep
            // every ion verbatim with plain base-loss semantics.
            return new ProductIonInterpretation(diagnosticIons, BaseLossType: null, BaseLossModification: null);
        }

        // The measured topology gives the primary ion its accurate monoisotopic m/z.
        var primaryIonIndex = ions.IndexOf(productIons.Max());
        if (primaryIonIndex >= 0)
        {
            ions[primaryIonIndex] = measured.Value.AccurateBaseCationMass;
        }

        return new ProductIonInterpretation(
            diagnosticIons,
            measured.Value.Behavior,
            measured.Value.Behavior == BaseLossBehavior.Default ? null : measured.Value.Formula);
    }

    /// <summary>
    /// The result of measuring how much of a modification sits on the base.
    /// </summary>
    /// <param name="Behavior">Default when nothing sits on the base; Modified otherwise.</param>
    /// <param name="Formula">The base-localized portion (null for Default).</param>
    /// <param name="AccurateBaseCationMass">The primary ion's accurate monoisotopic m/z: canonical base cation + portion.</param>
    private readonly record struct DerivedBaseLoss(BaseLossBehavior Behavior, ChemicalFormula? Formula, double AccurateBaseCationMass);

    /// <summary>
    /// Measures how much of the modification sits on the base, using the primary product ion: the
    /// largest published ion, which is the base cation (every other ion is that cation minus a small
    /// neutral such as water or ammonia, or a smaller sugar fragment).
    /// </summary>
    /// <param name="productIons">The published product ions.</param>
    /// <param name="canonicalBaseCationMass">[canonical free base + H]+, the unmodified residue's base cation.</param>
    /// <param name="modificationFormula">The modification's mass shift relative to the residue.</param>
    /// <returns>
    /// The uniquely measured portion: none of it (Default) for ribose-localized modifications such as
    /// 2'-O-methyls; a partial portion for split modifications such as N6,2'-O-dimethyladenosine
    /// (C1H2 of C2H4); the entire shift for base-localized modifications such as N6-methyladenosine.
    /// Null when the ions match no portion or two indistinguishable portions.
    /// </returns>
    private static DerivedBaseLoss? DeriveBaseLossModification(List<double> productIons, double canonicalBaseCationMass, ChemicalFormula modificationFormula)
    {
        var primaryIon = productIons.Max();

        // A portion is consistent with the measurement when adding it to the canonical base cation
        // reproduces the published ion at nominal precision.
        var matchingPortions = GetModificationPortions(modificationFormula)
            .Where(portion => SameNominalMz(canonicalBaseCationMass + portion.MonoisotopicMass, primaryIon))
            .ToList();

        if (matchingPortions.Count != 1)
        {
            // No portion matches, or two portions are isobaric at nominal precision: the ions do not
            // measure a unique answer, so keep the plain representation rather than guess.
            return null;
        }

        var baseLocalizedPortion = matchingPortions[0];
        var isEntirelySugarLocalized = string.IsNullOrEmpty(baseLocalizedPortion.Formula);
        return new DerivedBaseLoss(
            isEntirelySugarLocalized ? BaseLossBehavior.Default : BaseLossBehavior.Modified,
            isEntirelySugarLocalized ? null : baseLocalizedPortion,
            canonicalBaseCationMass + baseLocalizedPortion.MonoisotopicMass);
    }

    /// <summary>
    /// Every portion of the modification that could sit on the base: none of it, the entire shift, and —
    /// when the shift's element counts are all non-negative — each sub-formula in between.
    /// </summary>
    private static IEnumerable<ChemicalFormula> GetModificationPortions(ChemicalFormula modificationFormula)
    {
        // "None of the modification is on the base" is always a candidate.
        yield return new ChemicalFormula();

        var elements = ParseElementCounts(modificationFormula.Formula);
        var interiorIsEnumerable = elements.Count > 0
            && elements.All(e => e.Count >= 0)
            && elements.Aggregate(1L, (total, e) => total * (e.Count + 1)) <= 512;

        if (!interiorIsEnumerable)
        {
            // Shifts with negative counts (base conversions such as inosine's H-1N-1O1) cannot be split
            // into portions; for those, and for unusually large shifts, only the two extremes are
            // candidates, which already covers all-or-nothing topologies.
            yield return new ChemicalFormula(modificationFormula);
            yield break;
        }

        // For each element, a portion takes between zero and its full count in the modification.
        foreach (var portion in EnumeratePortions(elements, elementIndex: 0, prefix: string.Empty))
        {
            yield return portion;
        }
    }

    private static IEnumerable<ChemicalFormula> EnumeratePortions(List<(string Element, int Count)> elements, int elementIndex, string prefix)
    {
        if (elementIndex == elements.Count)
        {
            if (prefix.Length > 0)
            {
                yield return ChemicalFormula.ParseFormula(prefix);
            }

            yield break;
        }

        var (element, count) = elements[elementIndex];
        for (var taken = 0; taken <= count; taken++)
        {
            var withThisElement = taken == 0 ? prefix : prefix + element + taken;
            foreach (var portion in EnumeratePortions(elements, elementIndex + 1, withThisElement))
            {
                yield return portion;
            }
        }
    }

    private static List<(string Element, int Count)> ParseElementCounts(string formula)
    {
        // Hill notation: no spaces; counts of 1 are omitted; negative counts render as e.g. "H-1".
        var counts = new List<(string Element, int Count)>();
        foreach (Match match in Regex.Matches(formula, @"([A-Z][a-z]?)(-?\d+)?"))
        {
            var count = match.Groups[2].Success ? int.Parse(match.Groups[2].Value, CultureInfo.InvariantCulture) : 1;
            counts.Add((match.Groups[1].Value, count));
        }

        return counts;
    }

    private static bool SameNominalMz(double first, double second)
    {
        return (int)Math.Round(first, MidpointRounding.AwayFromZero) == (int)Math.Round(second, MidpointRounding.AwayFromZero);
    }

    private static Modification CreateModification(ModomicsDto dto, ModificationMotif motif, bool isTerminalCap,
        ChemicalFormula chemicalFormula, ProductIonInterpretation? productIonInterpretation)
    {
        var idString = dto.Id.ToString(CultureInfo.InvariantCulture);
        var modificationType = isTerminalCap ? "5' Terminal Cap" : "Modomics";
        var locationRestriction = isTerminalCap ? "5'-terminal." : "Anywhere.";
        var databaseReference = new Dictionary<string, IList<string>>
        {
            { "Modomics", new List<string> { dto.ShortName, dto.Name, idString } },
        };
        var keywords = new List<string> { dto.Abbrev, dto.ShortName, dto.Name };

        if (productIonInterpretation?.BaseLossType is { } baseLossType)
        {
            return new BaseModification(
                _originalId: dto.Name,
                _accession: idString,
                _modificationType: modificationType,
                _target: motif,
                _locationRestriction: locationRestriction,
                _chemicalFormula: chemicalFormula,
                _monoisotopicMass: chemicalFormula.MonoisotopicMass,
                _databaseReference: databaseReference,
                _keywords: keywords,
                _diagnosticIons: productIonInterpretation.DiagnosticIons,
                baseLossType: baseLossType,
                baseLossModification: productIonInterpretation.BaseLossModification);
        }

        return new Modification(
            _originalId: dto.Name,
            _accession: idString,
            _modificationType: modificationType,
            _target: motif,
            _locationRestriction: locationRestriction,
            _chemicalFormula: chemicalFormula,
            _monoisotopicMass: chemicalFormula.MonoisotopicMass,
            _databaseReference: databaseReference,
            _keywords: keywords,
            _diagnosticIons: productIonInterpretation?.DiagnosticIons);
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
