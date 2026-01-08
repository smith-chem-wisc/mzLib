using Chemistry;
using Omics;
using Omics.Modifications;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Reflection;
using System.Security.Policy;
using UsefulProteomicsDatabases;
using static TorchSharp.torch.optim.lr_scheduler.impl.CyclicLR;

namespace Readers;

public static class ModificationConverter
{
    #region All Known Mods for MetaMorpheus - Use this to check for known mods to compare against 

    // Cache of previously successful modification conversions to reduce the number of times the same modification is converted
    private static readonly ConcurrentDictionary<(string, char), Modification> ModificationCache;

    public static List<Modification> AllKnownMods { get; }
    public static Dictionary<string, Modification> AllModsKnown { get; }

    // Separate lists for different mod types
    private static List<Modification> _uniprotMods;
    private static List<Modification> _metamorpheusMods;

    static ModificationConverter()
    {
        ModificationCache = new ConcurrentDictionary<(string, char), Modification>();

        var info = Assembly.GetExecutingAssembly().GetName();
        var name = info.Name;

        var unimodStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.unimod.xml");
        var unimodMods = UnimodLoader.ReadMods(unimodStream);

        var psiModStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.PSI-MOD.obo.xml");
        var psiModObo = Loaders.LoadPsiMod(psiModStream);
        var formalChargeDict = Loaders.GetFormalChargesDictionary(psiModObo);

        var uniprotStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.ptmlist.txt");
        var uniprotPtms = PtmListLoader.ReadModsFromFile(new StreamReader(uniprotStream), formalChargeDict, out _).ToList();

        var modsTextStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.Mods.txt");
        var modsTextMods = PtmListLoader.ReadModsFromFile(new StreamReader(modsTextStream), formalChargeDict, out _).ToList();

        // Store UniProt mods separately (they have ModificationType == "UniProt")
        _uniprotMods = uniprotPtms.Where(m => m.ModificationType == "UniProt").ToList();

        // MetaMorpheus mods are everything else
        _metamorpheusMods = modsTextMods;

        AllKnownMods = unimodMods.Concat(uniprotPtms).Concat(modsTextMods).ToList();
        AllModsKnown = AllKnownMods
            .DistinctBy(m => m.IdWithMotif)
            .ToDictionary(m => m.IdWithMotif);
    }

    #endregion

    public static Dictionary<int, Modification> GetModificationDictionaryFromFullSequence(string fullSequence, Dictionary<string, Modification>? allModsKnown = null)
    {
        var allModsOneIsNterminus = new Dictionary<int, Modification>();
        allModsKnown ??= AllModsKnown;
        var baseSequence =  IBioPolymerWithSetMods.GetBaseSequenceFromFullSequence(fullSequence);
        int currentModStart = 0;
        int currentModificationLocation = 1;
        bool currentlyReadingMod = false;
        int bracketCount = 0;

        for (int r = 0; r < fullSequence.Length; r++)
        {
            char c = fullSequence[r];
            if (c == '[')
            {
                currentlyReadingMod = true;
                if (bracketCount == 0)
                {
                    currentModStart = r + 1;
                }
                bracketCount++;
            }
            else if (c == ']')
            {
                string modId = null;
                bracketCount--;
                if (bracketCount == 0)
                {
                    Modification? mod = null;
                    Exception? toThrow = null;
                    string modString = fullSequence.Substring(currentModStart, r - currentModStart);
                    try
                    {
                        //remove the beginning section (e.g. "Fixed", "Variable", "Uniprot") if present
                        int splitIndex = modString.IndexOf(':');
                        modId = splitIndex > 0 ? modString.Substring(splitIndex + 1, modString.Length - splitIndex - 1) : modString;

                        if (!allModsKnown.TryGetValue(modId, out mod))
                        {
                            toThrow = new MzLibUtil.MzLibException(
                                "Could not find modification while reading string: " + fullSequence);
                        }
                    }
                    catch (Exception e)
                    {
                        toThrow = new MzLibUtil.MzLibException(
                            "Error while trying to parse string into peptide: " + e.Message, e);
                    }

                    // Not standard MM format
                    if (mod == null && toThrow != null)
                    {
                        try
                        {
                            char modifiedResidue;
                            if (currentModificationLocation - 2 < 0) // N-Terminal
                                modifiedResidue = 'X';
                            else
                                modifiedResidue = baseSequence[currentModificationLocation - 2];
                            mod = GetClosestMod(modString, modifiedResidue);
                        }
                        catch (Exception e)
                        {
                            toThrow = new MzLibUtil.MzLibException(
                                "Could not find modification while reading string: " + fullSequence +
                                Environment.NewLine + "Also could not find closest match for modification: " + modId +
                                Environment.NewLine + e.Message, e);
                        }
                    }

                    if (mod == null)
                    {
                        throw toThrow!;
                    }

                    // Set the C-terminus modification index to its OneIsNTerminus Index.
                    // Checks if the location restriction for the mod contains C-terminal' (for protein and peptide BioPolymer objects)
                    // or '3'-terminal' (for nucleic acid BioPolymer objects) and if we are at the last residue of the full sequence.
                    if ((mod.LocationRestriction.Contains("C-terminal.") || mod.LocationRestriction.Contains("3'-terminal.") && r == fullSequence.Length - 1))
                    {
                        currentModificationLocation = baseSequence.Length + 2;
                    }

                    allModsOneIsNterminus.Add(currentModificationLocation, mod);
                    currentlyReadingMod = false;
                }
            }
            else if (!currentlyReadingMod && c != '-')
            {
                currentModificationLocation++;
            }
            //else do nothing
        }

        return allModsOneIsNterminus;
    }


    /// <summary>
    /// Gets the closest modification from the list of all known modifications that matches the given localized modification,
    /// with a preference for the specified naming convention
    /// </summary>
    public static Modification? GetClosestMod(Modification mod, char modifiedResidue, ModificationNamingConvention convention = ModificationNamingConvention.MetaMorpheus, IHasChemicalFormula? chemicalFormula = null)
    {
        var modsToSearch = GetModsByNamingConvention(convention);

        if (modsToSearch.Contains(mod))
            return mod;

        try
        {
            return GetClosestMod(mod.OriginalId, modifiedResidue, modsToSearch, chemicalFormula);
        }
        catch (KeyNotFoundException e) 
        { 
            
        }

        return null;
    }

    /// <summary>
    /// Gets the closest modification from the list of all known modifications that matches the given localized modification.
    /// </summary>
    public static Modification GetClosestMod(string name, char modifiedResidue, IList<Modification>? allKnownMods = null, IHasChemicalFormula? chemicalFormula = null)
    {
        allKnownMods ??= AllKnownMods;
        var cacheKey = (name, modifiedResidue);
        // if we have done this conversion before, just return it. 
        if (ModificationCache.TryGetValue(cacheKey, out var cachedModification))
        {
            return cachedModification;
        }

        if (chemicalFormula != null)
            allKnownMods = allKnownMods.Where(mod => mod.ChemicalFormula.Equals(chemicalFormula.ThisChemicalFormula)).ToList();

        // get rid of this annoying -L- suffix that is added to some mods in UniProt
        var trimmedName = name.Split(new[] { "-L-" }, StringSplitOptions.None)[0];
        var nameContaining = allKnownMods.Where(p => p.IdWithMotif.Contains(trimmedName)).ToHashSet();
        if (nameContaining.Count == 0 && trimmedName.Contains("ation"))
            nameContaining = allKnownMods.Where(p => p.IdWithMotif.Contains(trimmedName.Replace("ation", ""))).ToHashSet();

        var motifMatching = allKnownMods.Where(p => p.IdWithMotif.Contains($" on {modifiedResidue}") || p.IdWithMotif.Contains(" on X")).ToHashSet();

        var goodMatches = nameContaining.Intersect(motifMatching).ToList();
        switch (goodMatches.Count)
        {
            // if exact match by name with no ambiguity, return it
            case 1:
                cachedModification = goodMatches.First();
                break;

            // Many matched by name and motif, see if we can whittle down possible options
            case > 1:
                // remove those that are labels if we have other options
                if (goodMatches.Count(p => p.IdWithMotif.StartsWith("Label")) < goodMatches.Count)
                    foreach (var mod in goodMatches.Where(p => p.IdWithMotif.StartsWith("Label")).ToList())
                        goodMatches.Remove(mod);

                // TODO: anything else that we can do to reduce options here? 
                
                if (goodMatches.Count == 1) // if we have only one left, return it
                    cachedModification = goodMatches.First();
                else // if we still have many, we need to do some scoring
                    cachedModification = goodMatches
                        .OrderByDescending(mod => GetOverlapScore(mod.IdWithMotif, trimmedName))
                        .ThenByDescending(p => p.Target.ToString() == modifiedResidue.ToString())
                        .ThenByDescending(p => p.ModificationType.Length)
                        .ThenBy(p => p.IdWithMotif.Length).First();
                break;

            // None matched by name and motif, Calculate overlap score by substring overlap of al possible matches. 
            // and select the modification with the highest score, better matching residue, then shortest mod name in case of tie
            case < 1:
                var candidates = nameContaining.Union(motifMatching).ToList();
                if (!candidates.Any())
                {
                    throw new KeyNotFoundException($"Could not find a modification that matches {name} on {modifiedResidue}");
                }
                
                cachedModification = candidates
                    .OrderByDescending(mod => GetOverlapScore(mod.IdWithMotif, trimmedName))
                    .ThenByDescending(p => p.Target.ToString() == modifiedResidue.ToString())
                    .ThenByDescending(p => p.ModificationType.Length)
                    .ThenBy(p => p.IdWithMotif.Length).First();
                break;
        }

        if (cachedModification == null)
            throw new KeyNotFoundException($"Could not find a modification that matches {name} on {modifiedResidue}");

        ModificationCache[cacheKey] = cachedModification;
        return cachedModification;
    }

    /// <summary>
    /// Calculates the overlap score between the modification ID with motif and the trimmed name.
    /// The score represents the length of the longest common substring between the two strings.
    /// </summary>
    /// <param name="idWithMotif">The modification ID with motif.</param>
    /// <param name="trimmedName">The trimmed name of the modification.</param>
    /// <returns>The overlap score, which is the length of the longest common substring.</returns>
    private static int GetOverlapScore(string idWithMotif, string trimmedName)
    {
        int overlapScore = 0;
        for (int i = 0; i < idWithMotif.Length; i++)
        {
            for (int j = 0; j < trimmedName.Length; j++)
            {
                int k = 0;
                while (i + k < idWithMotif.Length && j + k < trimmedName.Length && idWithMotif[i + k] == trimmedName[j + k])
                {
                    k++;
                }
                overlapScore = Math.Max(overlapScore, k);
            }
        }
        return overlapScore;
    }

    private static List<Modification> GetModsByNamingConvention(ModificationNamingConvention convention)
    {
        return convention switch
        {
            ModificationNamingConvention.MetaMorpheus => _metamorpheusMods,
            ModificationNamingConvention.UniProt => _uniprotMods,
            ModificationNamingConvention.Mixed => AllKnownMods,
            _ => AllKnownMods,
        };
    }

    public static void ConvertMods(this IBioPolymerWithSetMods withSetMods, ModificationNamingConvention targetConvention)
    {
        foreach (var kvp in withSetMods.AllModsOneIsNterminus)
        {
            char residue;

            // N-Terminal
            if (kvp.Key == 1)
                residue = withSetMods.BaseSequence[0];
            // C-Terminal
            else if (kvp.Key == withSetMods.Length + 2)
                residue = withSetMods.BaseSequence[withSetMods.Length - 1];
            else
                residue = withSetMods.BaseSequence[kvp.Key - 2];

            var convertedMod = GetClosestMod(kvp.Value, residue, targetConvention);

            if (convertedMod is null)
                withSetMods.AllModsOneIsNterminus.Remove(kvp.Key);
            else if (convertedMod.ChemicalFormula.Equals(kvp.Value.ChemicalFormula))
                withSetMods.AllModsOneIsNterminus[kvp.Key] = convertedMod;
            else
                Debugger.Break();
        }
    }

    public static void ConvertMods(this IBioPolymer bioPolymer, ModificationNamingConvention targetConvention)
    {
        foreach (var kvp in bioPolymer.OneBasedPossibleLocalizedModifications)
        {
            char residue;

            // N-Terminal
            if (kvp.Key == 1)
                residue = bioPolymer.BaseSequence[0];
            // C-Terminal
            else if (kvp.Key == bioPolymer.Length + 2)
                residue = bioPolymer.BaseSequence[bioPolymer.Length - 1];
            else
                residue = bioPolymer.BaseSequence[kvp.Key - 1];

            for (var index = 0; index < kvp.Value.Count; index++)
            {
                var mod = kvp.Value[index];

                var convertedMod = GetClosestMod(mod, residue, targetConvention, mod.ChemicalFormula);

                if (convertedMod is null)
                {
                    bioPolymer.OneBasedPossibleLocalizedModifications[kvp.Key].RemoveAt(index);
                    if (bioPolymer.OneBasedPossibleLocalizedModifications[kvp.Key].Count == 0)
                        bioPolymer.OneBasedPossibleLocalizedModifications.Remove(kvp.Key);
                }
                else if (convertedMod.Equals(mod) && convertedMod.FeatureType != "UniProt")
                    Debugger.Break();
                else if (convertedMod.ChemicalFormula.Equals(mod.ChemicalFormula))
                    kvp.Value[index] = convertedMod;
                else
                    Debugger.Break();
            }
        }

        foreach (var kvp in bioPolymer.OriginalNonVariantModifications)
        {
            char residue;

            // N-Terminal
            if (kvp.Key == 1)
                residue = bioPolymer.BaseSequence[0];
            // C-Terminal
            else if (kvp.Key == bioPolymer.Length + 2)
                residue = bioPolymer.BaseSequence[bioPolymer.Length - 1];
            else
                residue = bioPolymer.BaseSequence[kvp.Key - 1];

            for (var index = 0; index < kvp.Value.Count; index++)
            {
                var mod = kvp.Value[index];
                var convertedMod = GetClosestMod(mod, residue, targetConvention, mod.ChemicalFormula);

                if (convertedMod is null)
                {
                    bioPolymer.OriginalNonVariantModifications[kvp.Key].RemoveAt(index);
                    if (bioPolymer.OriginalNonVariantModifications[kvp.Key].Count == 0)
                        bioPolymer.OriginalNonVariantModifications.Remove(kvp.Key);
                }
                else if (convertedMod.ChemicalFormula.Equals(mod.ChemicalFormula))
                    kvp.Value[index] = convertedMod;
                else
                    Debugger.Break();
            }
        }

        // Convert mods in SequenceVariations
        // Note: SequenceVariation.OneBasedModifications are indexed relative to the variant sequence, not the protein
        foreach (var variant in bioPolymer.SequenceVariations)
        {
            foreach (var kvp in variant.OneBasedModifications)
            {
                for (var index = 0; index < kvp.Value.Count; index++)
                {
                    var mod = kvp.Value[index];
                    char residue;

                    // Determine residue from the variant sequence
                    // Keys in variant.OneBasedModifications are relative to the variant sequence
                    if (string.IsNullOrEmpty(variant.VariantSequence))
                    {
                        // If no variant sequence, use 'X' as unknown
                        residue = 'X';
                    }
                    else if (kvp.Key == 1 && variant.VariantSequence.Length > 0)
                    {
                        // First position in variant
                        residue = variant.VariantSequence[0];
                    }
                    else if (kvp.Key <= variant.VariantSequence.Length)
                    {
                        // Internal position in variant (convert from 1-based to 0-based)
                        residue = variant.VariantSequence[kvp.Key - 1];
                    }
                    else
                    {
                        // Position beyond variant sequence
                        residue = 'X';
                    }

                    var convertedMod = GetClosestMod(mod, residue, targetConvention, mod.ChemicalFormula);
                    if (convertedMod is null)
                        kvp.Value.RemoveAt(index);
                    else if (convertedMod.ChemicalFormula.Equals(mod.ChemicalFormula))
                        kvp.Value[index] = convertedMod;
                    else
                        Debugger.Break();
                }
            }
        }

        // Convert mods in AppliedSequenceVariations
        foreach (var variant in bioPolymer.AppliedSequenceVariations)
        {
            foreach (var kvp in variant.OneBasedModifications)
            {
                for (var index = 0; index < kvp.Value.Count; index++)
                {
                    var mod = kvp.Value[index];
                    char residue;

                    // Determine residue from the variant sequence
                    // Keys in variant.OneBasedModifications are relative to the variant sequence
                    if (string.IsNullOrEmpty(variant.VariantSequence))
                    {
                        // If no variant sequence, use 'X' as unknown
                        residue = 'X';
                    }
                    else if (kvp.Key == 1 && variant.VariantSequence.Length > 0)
                    {
                        // First position in variant
                        residue = variant.VariantSequence[0];
                    }
                    else if (kvp.Key <= variant.VariantSequence.Length)
                    {
                        // Internal position in variant (convert from 1-based to 0-based)
                        residue = variant.VariantSequence[kvp.Key - 1];
                    }
                    else
                    {
                        // Position beyond variant sequence
                        residue = 'X';
                    }

                    var convertedMod = GetClosestMod(mod, residue, targetConvention, mod.ChemicalFormula);
                    if (convertedMod is null)
                        kvp.Value.RemoveAt(index);
                    else if (convertedMod.ChemicalFormula.Equals(mod.ChemicalFormula))
                        kvp.Value[index] = convertedMod;
                    else
                        Debugger.Break();
                }
            }
        }
    }
}