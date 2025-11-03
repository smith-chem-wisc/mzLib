using Omics.Modifications;
using System.Collections.Concurrent;
using System.Reflection;
using Omics;
using UsefulProteomicsDatabases;

namespace Readers;

public static class ModificationConverter
{
    #region All Known Mods for MetaMorpheus - Use this to check for known mods to compare against 

    // Cache of previously successful modification conversions to reduce the number of times the same modification is converted
    private static readonly ConcurrentDictionary<(string, char), Modification> _modificationCache;

    internal static List<Modification> AllKnownMods;
    internal static Dictionary<string, Modification> AllModsKnown;

    static ModificationConverter()
    {
        _modificationCache = new ConcurrentDictionary<(string, char), Modification>();

        var info = Assembly.GetExecutingAssembly().GetName();
        var name = info.Name;

        var unimodStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.unimod.xml");
        var unimodMods = UnimodLoader.ReadMods(unimodStream);

        var psiModStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.PSI-MOD.obo.xml");
        var psiModObo = Loaders.LoadPsiMod(psiModStream);
        var formalChargeDict = Loaders.GetFormalChargesDictionary(psiModObo);

        var uniprotStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.ptmlist.txt");
        var uniprotMods = PtmListLoader.ReadModsFromFile(new StreamReader(uniprotStream), formalChargeDict, out _);

        var modsTextStream = Assembly.GetExecutingAssembly().GetManifestResourceStream($"{name}.Resources.Mods.txt");
        var modsTextMods = PtmListLoader.ReadModsFromFile(new StreamReader(modsTextStream), formalChargeDict, out _);

        AllKnownMods = unimodMods.Concat(uniprotMods).Concat(modsTextMods).ToList();
        AllModsKnown = AllKnownMods
            .DistinctBy(m => (m.ModificationType, m.IdWithMotif))
            .ToDictionary(m => $"{m.ModificationType}:{ m.IdWithMotif }");
    }

    #endregion


    public static Dictionary<int, Modification> GetModificationDictionaryFromFullSequence(string fullSequence,
    Dictionary<string, Modification>? allModsKnown = null)
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
                        //remove the beginning section (e.g. "Fixed", "Variable", "Uniprot")
                        int splitIndex = modString.IndexOf(':');
                        string modType = modString.Substring(0, splitIndex);
                        modId = modString.Substring(splitIndex + 1, modString.Length - splitIndex - 1);

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
    /// Gets the closest modification from the list of all known modifications that matches the given localized modification.
    /// </summary>
    public static Modification GetClosestMod(string name, char modifiedResidue, IList<Modification>? allKnownMods = null)
    {
        allKnownMods ??= AllKnownMods;
        var cacheKey = (name, modifiedResidue);
        // if we have done this conversion before, just return it. 
        if (_modificationCache.TryGetValue(cacheKey, out var cachedModification))
        {
            return cachedModification;
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
}