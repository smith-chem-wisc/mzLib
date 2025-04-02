using Omics.Modifications;
using System.Collections.Concurrent;
using System.Reflection;
using UsefulProteomicsDatabases;

namespace Readers;

public interface ILocalizedModification
{
    public string Name { get; }
    public int OneBasedLocalization { get; }
    public char ModifiedResidue { get; }
    public double MonoisotopicMass { get; }
}

public static class ModificationExtensions
{
    #region All Known Mods for MetaMorpheus - Use this to check for known mods to compare against 

    // Cache of previously successful modification conversions to reduce the number of times the same modification is converted
    private static readonly ConcurrentDictionary<(string, char), Modification> _modificationCache;

    internal static List<Modification> AllKnownMods;

    static ModificationExtensions()
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
    }

    #endregion

    /// <summary>
    /// Get the string representation of the modification in MetaMorpheus format
    /// </summary>
    public static string GetMetaMorpheusFullSequenceString(this ILocalizedModification mod, IList<Modification>? allKnownMods = null)
    {
        var mmMod = mod.GetClosestMod(allKnownMods);

        string category = mmMod.ModificationType switch
        {
            "Unimod" when mmMod.OriginalId.Contains("Carbamido") => "Common Fixed",
            "Unimod" when mmMod.OriginalId.Contains("Oxidation") => "Common Variable",
            "Unimod" when mmMod.OriginalId.Contains("Phosphoryl") => "Common Biological",
            "Unimod" when mmMod.OriginalId.Contains("Acetyl") => "Common Biological",
            "Unimod" when mmMod.OriginalId.Contains("Methy") => "Common Biological",
            _ => mmMod.ModificationType
        };

        return $"[{category}:{mmMod.OriginalId} on {mmMod.Target}]";
    }

    /// <summary>
    /// Gets the closest modification from the list of all known modifications that matches the given localized modification.
    /// </summary>
    public static Modification GetClosestMod(this ILocalizedModification modToMatch, IList<Modification>? allKnownMods = null)
    {
        var name = modToMatch.Name;
        var modifiedResidue = modToMatch.ModifiedResidue;
        return GetClosestMod(name, modifiedResidue, allKnownMods);
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

        // Gather those that match by name and Motif
        var matching = allKnownMods.Where(p =>
                p.IdWithMotif.Contains(trimmedName)
                && (p.IdWithMotif.Contains($" on {modifiedResidue}") || p.IdWithMotif.Contains(" on X")))
            .ToList();

        switch (matching.Count)
        {
            // if exact match by name with no ambiguity, return it
            case 1:
                cachedModification = matching[0];
                break;
            // if none matched by name alone, add all that have the desired motif to matching set
            case < 1:
                {
                    var motifMatching = allKnownMods.Where(p =>
                        p.IdWithMotif.Contains($" on {modifiedResidue}") || p.IdWithMotif.Contains(" on X"));
                    matching.AddRange(motifMatching);
                    break;
                }
        }

        // if multiple match by name and motif, but all have the same name, return the one with the correct motif 
        if (matching.Count > 1 && matching.DistinctBy(p => p.OriginalId).Count() == 1)
        {
            var exactMatch = matching.FirstOrDefault(p => p.IdWithMotif.Contains($" on {modifiedResidue}"));
            if (exactMatch is not null)
                cachedModification = exactMatch;
            else
            {
                var ambiguousMatch = matching.FirstOrDefault(p => p.IdWithMotif.Contains(" on X"));
                if (ambiguousMatch is not null)
                    cachedModification = ambiguousMatch;
            }
        }

        // if nothing above worked, Calculate overlap score by substring overlap and select the modification with the highest score
        if (cachedModification is null)
        {
            // Calculate overlap score and select the modification with the highest score
            cachedModification = matching.MaxBy(mod => GetOverlapScore(mod.IdWithMotif, trimmedName))!;
        }

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