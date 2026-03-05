using Omics;
using Omics.Modifications;
using Omics.Modifications.Conversion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace Readers;

public static class ModificationConverter
{
    #region All Known Mods for MetaMorpheus - Use this to check for known mods to compare against 

    private static readonly ConcurrentDictionary<(string, char), Modification> _modificationCache;
    private static readonly ConcurrentDictionary<(string key, char residue), Modification> _explicitIdMatchCache;
    private static readonly object _cacheLock = new();
    private static int _cacheVersion = -1;
    private static ModificationCrossRefIndex _crossRefIndex = ModificationCrossRefIndex.Global;

    internal static List<Modification> AllKnownMods => Mods.AllProteinModsList;
    internal static Dictionary<string, Modification> AllModsKnown => Mods.AllKnownProteinModsDictionary;

    private static readonly Regex UnimodRegex = new("(?i)UNIMOD:\\s*(\\d+)", RegexOptions.Compiled);
    private static readonly Regex PsiModFullRegex = new("(?i)PSI-MOD:\\s*(MOD:\\d+)", RegexOptions.Compiled);
    private static readonly Regex PsiModShortRegex = new("(?i)\bMOD:\\d+\b", RegexOptions.Compiled);
    private static readonly Regex ResidRegex = new("(?i)RESID:\\s*([A-Z0-9]+)", RegexOptions.Compiled);
    private static readonly Regex PtmRegex = new("(?i)PTM-\\d+", RegexOptions.Compiled);

    static ModificationConverter()
    {
        _modificationCache = new ConcurrentDictionary<(string, char), Modification>();
        _explicitIdMatchCache = new ConcurrentDictionary<(string, char), Modification>();
        _cacheVersion = Mods.RegistryVersion;
    }

    #endregion


    public static Dictionary<int, Modification> GetModificationDictionaryFromFullSequence(string fullSequence,
    Dictionary<string, Modification>? allModsKnown = null)
    {
        var allModsOneIsNterminus = new Dictionary<int, Modification>();
        allModsKnown ??= AllModsKnown;
        EnsureCachesCurrent();
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
                        char modifiedResidue = currentModificationLocation - 2 < 0
                            ? 'X'
                            : baseSequence[currentModificationLocation - 2];

                        if (TryResolveByExplicitIds(modifiedResidue, out var resolvedFromId, modString, modId))
                        {
                            mod = resolvedFromId;
                            toThrow = null;
                        }
                    }

                    // Still not resolved, fall back to heuristic search
                    if (mod == null && toThrow != null)
                    {
                        char modifiedResidue;
                        if (currentModificationLocation - 2 < 0) // N-Terminal
                            modifiedResidue = 'X';
                        else
                            modifiedResidue = baseSequence[currentModificationLocation - 2];

                        try
                        {
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
    /// Gets the closest modification from the list of all known modifications that matches the given localized modification.
    /// </summary>
    public static Modification GetClosestMod(string name, char modifiedResidue, IList<Modification>? allKnownMods = null)
    {
        EnsureCachesCurrent();
        allKnownMods ??= AllKnownMods;
        var cacheKey = (name, modifiedResidue);
        // if we have done this conversion before, just return it. 
        if (_modificationCache.TryGetValue(cacheKey, out var cachedModification))
        {
            return cachedModification;
        }

        if (TryResolveByExplicitIds(modifiedResidue, out var resolvedById, name))
        {
            _modificationCache[cacheKey] = resolvedById;
            return resolvedById;
        }

        // get rid of this annoying -L- suffix that is added to some mods in UniProt
        var trimmedName = name.Split(new[] { "-L-" }, StringSplitOptions.None)[0];

        var nameContaining = allKnownMods.Where(p => p.IdWithMotif.Contains(trimmedName)).ToHashSet();
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
                cachedModification = nameContaining.Union(motifMatching)
                    .OrderByDescending(mod => GetOverlapScore(mod.IdWithMotif, trimmedName))
                    .ThenByDescending(p => p.Target.ToString() == modifiedResidue.ToString())
                    .ThenByDescending(p => p.ModificationType.Length)
                    .ThenBy(p => p.IdWithMotif.Length).First();
                break;
        }

        if (cachedModification == null)
            throw new KeyNotFoundException($"Could not find a modification that matches {name} on {modifiedResidue}");

        _modificationCache[cacheKey] = cachedModification;
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

    private static bool TryResolveByExplicitIds(char modifiedResidue, out Modification match, params string?[] annotations)
    {
        match = null;
        if (annotations == null)
        {
            return false;
        }

        foreach (var annotation in annotations)
        {
            if (string.IsNullOrWhiteSpace(annotation))
            {
                continue;
            }

            foreach (var reference in ExtractDatabaseReferences(annotation))
            {
                var cacheKey = (reference.cacheKey, modifiedResidue);
                if (_explicitIdMatchCache.TryGetValue(cacheKey, out match))
                {
                    return true;
                }

                if (TryResolveReference(reference.database, reference.accession, modifiedResidue, out match))
                {
                    _explicitIdMatchCache[cacheKey] = match;
                    return true;
                }
            }
        }

        return false;
    }

    private static IEnumerable<(string? database, string accession, string cacheKey)> ExtractDatabaseReferences(string annotation)
    {
        var seen = new HashSet<string>(StringComparer.OrdinalIgnoreCase);

        foreach (Match match in UnimodRegex.Matches(annotation))
        {
            var accession = match.Groups[1].Value.Trim();
            var cacheKey = $"UNIMOD:{accession.ToUpperInvariant()}";
            if (seen.Add(cacheKey))
            {
                yield return ("UNIMOD", accession, cacheKey);
            }
        }

        foreach (Match match in PsiModFullRegex.Matches(annotation))
        {
            var accession = match.Groups[1].Value.Trim();
            var cacheKey = $"PSI-MOD:{accession.ToUpperInvariant()}";
            if (seen.Add(cacheKey))
            {
                yield return ("PSI-MOD", accession, cacheKey);
            }
        }

        foreach (Match match in PsiModShortRegex.Matches(annotation))
        {
            var accession = match.Value.Trim();
            var cacheKey = $"PSI-MOD:{accession.ToUpperInvariant()}";
            if (seen.Add(cacheKey))
            {
                yield return ("PSI-MOD", accession, cacheKey);
            }
        }

        foreach (Match match in ResidRegex.Matches(annotation))
        {
            var accession = match.Groups[1].Value.Trim();
            var cacheKey = $"RESID:{accession.ToUpperInvariant()}";
            if (seen.Add(cacheKey))
            {
                yield return ("RESID", accession, cacheKey);
            }
        }

        foreach (Match match in PtmRegex.Matches(annotation))
        {
            var accession = match.Value.Trim().ToUpperInvariant();
            var cacheKey = $"ACC:{accession}";
            if (seen.Add(cacheKey))
            {
                yield return (null, accession, cacheKey);
            }
        }
    }

    private static bool TryResolveReference(string? databaseName, string accession, char modifiedResidue, out Modification match)
    {
        IReadOnlyList<Modification> candidates = databaseName == null
            ? _crossRefIndex.GetByAccession(accession)
            : _crossRefIndex.GetByDatabaseId(databaseName, accession);

        if (candidates.Count == 0)
        {
            match = null!;
            return false;
        }

        if (candidates.Count == 1)
        {
            match = candidates[0];
            return true;
        }

        var residueMatches = candidates.Where(c => ResidueMatches(c, modifiedResidue)).ToList();
        if (residueMatches.Count == 1)
        {
            match = residueMatches[0];
            return true;
        }

        match = null!;
        return false;
    }

    private static bool ResidueMatches(Modification mod, char residue)
    {
        if (residue == 'X' || mod.Target == null)
        {
            return true;
        }

        var motif = mod.Target.ToString();
        if (string.IsNullOrEmpty(motif))
        {
            return true;
        }

        var targetResidue = char.ToUpperInvariant(residue);
        foreach (var ch in motif)
        {
            if (!char.IsLetter(ch))
            {
                continue;
            }

            var upper = char.ToUpperInvariant(ch);
            if (upper == 'X' || upper == targetResidue)
            {
                return true;
            }
        }

        return false;
    }

    private static void EnsureCachesCurrent()
    {
        var currentVersion = Mods.RegistryVersion;
        if (_cacheVersion == currentVersion)
        {
            return;
        }

        lock (_cacheLock)
        {
            if (_cacheVersion == currentVersion)
            {
                return;
            }

            _modificationCache.Clear();
            _explicitIdMatchCache.Clear();
            _crossRefIndex = ModificationCrossRefIndex.Global;
            _cacheVersion = currentVersion;
        }
    }
}
