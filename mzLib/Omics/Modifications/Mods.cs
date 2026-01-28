using Omics.Modifications.IO;
using System.Reflection;

namespace Omics.Modifications;

public static class Mods
{
    private static readonly Lazy<Dictionary<string, Modification>> _allKnownProteinMods;
    private static readonly Lazy<List<Modification>> _allProteinModsList;
    private static readonly Lazy<Dictionary<string, Modification>> _allKnownRnaMods;
    private static readonly Lazy<List<Modification>> _allRnaModsList;

    private static readonly object _cacheLock = new object();

    static Mods()
    {
        _allKnownProteinMods = new Lazy<Dictionary<string, Modification>>(LoadAllProteinModifications);
        _allProteinModsList = new Lazy<List<Modification>>(() => _allKnownProteinMods.Value.Values.ToList());
        _allKnownRnaMods = new Lazy<Dictionary<string, Modification>>(LoadAllRnaModifications);
        _allRnaModsList = new Lazy<List<Modification>>(() => _allKnownRnaMods.Value.Values.ToList());
    }

    #region Public Properties

    /// <summary>
    /// All known protein modifications indexed by IdWithMotif
    /// </summary>
    public static Dictionary<string, Modification> AllKnownProteinMods => _allKnownProteinMods.Value;

    /// <summary>
    /// All known protein modifications as a list
    /// </summary>
    public static List<Modification> AllProteinModsList => _allProteinModsList.Value;

    /// <summary>
    /// All known RNA modifications indexed by IdWithMotif
    /// </summary>
    public static Dictionary<string, Modification> AllKnownRnaMods => _allKnownRnaMods.Value;

    /// <summary>
    /// All known RNA modifications as a list
    /// </summary>
    public static List<Modification> AllRnaModsList => _allRnaModsList.Value;

    /// <summary>
    /// Combined list of all known modifications (protein + RNA)
    /// </summary>
    public static List<Modification> AllKnownMods =>
        AllProteinModsList.Concat(AllRnaModsList).ToList();

    /// <summary>
    /// Combined dictionary of all known modifications (protein + RNA)
    /// RNA mods take precedence in case of conflicts
    /// </summary>
    public static Dictionary<string, Modification> AllModsKnown
    {
        get
        {
            var combined = new Dictionary<string, Modification>(AllKnownProteinMods);
            foreach (var kvp in AllKnownRnaMods)
            {
                combined[kvp.Key] = kvp.Value;
            }
            return combined;
        }
    }

    #endregion

    #region Loading Methods

    /// <summary>
    /// Loads Protein modifications from embedded resources and custom sources
    /// </summary>
    private static Dictionary<string, Modification> LoadAllProteinModifications()
    {
        var allMods = new Dictionary<string, Modification>();
        var assembly = Assembly.GetExecutingAssembly();
        var assemblyName = assembly.GetName().Name;

        try
        {
            // 1. Load Unimod
            var unimodStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.unimod.xml");
            if (unimodStream != null)
            {
                var unimodMods = ModificationLoader.ReadModsFromUnimod(unimodStream);
                AddModsToDictionary(allMods, unimodMods, "Unimod");
            }

            // 2. Load PSI-MOD and get formal charges
            Dictionary<string, int> formalChargeDict = new Dictionary<string, int>();
            var psiModStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.PSI-MOD.obo.xml");
            if (psiModStream != null)
            {
                var psiModObo = ModificationLoader.LoadPsiMod(psiModStream);
                formalChargeDict = ModificationLoader.GetFormalChargesDictionary(psiModObo);
            }

            // 3. Load UniProt ptmlist
            var uniprotStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.ptmlist.txt");
            if (uniprotStream != null)
            {
                using (var reader = new StreamReader(uniprotStream))
                {
                    var uniprotMods = ModificationLoader.ReadModsFromFile(reader, formalChargeDict,
                        out var filteredWarnings);
                    AddModsToDictionary(allMods, uniprotMods, "UniProt");

                    // Optionally log filtered mods
                    foreach (var (mod, warning) in filteredWarnings)
                    {
                        System.Diagnostics.Debug.WriteLine($"Filtered UniProt mod: {warning}");
                    }
                }
            }

            // 4. Load custom Mods.txt
            var modsTextStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.Mods.txt");
            if (modsTextStream != null)
            {
                using (var reader = new StreamReader(modsTextStream))
                {
                    var modsTextMods = ModificationLoader.ReadModsFromFile(reader, formalChargeDict,
                        out var filteredWarnings);
                    AddModsToDictionary(allMods, modsTextMods, "Custom");

                    foreach (var (mod, warning) in filteredWarnings)
                    {
                        System.Diagnostics.Debug.WriteLine($"Filtered custom mod: {warning}");
                    }
                }
            }
        }
        catch (Exception ex)
        {
            System.Diagnostics.Debug.WriteLine($"Error loading protein modifications: {ex.Message}");
            throw new InvalidOperationException("Failed to load protein modifications from embedded resources", ex);
        }

        return allMods;
    }


    /// <summary>
    /// Loads RNA modifications from embedded resources and custom sources
    /// </summary>
    private static Dictionary<string, Modification> LoadAllRnaModifications()
    {
        var allMods = new Dictionary<string, Modification>();
        var assembly = Assembly.GetExecutingAssembly();
        var assemblyName = assembly.GetName().Name;

        try
        {
            // Load RNA-specific mods if we have a resource for them
            // For now, this is a placeholder for future RNA mod sources
            var rnaModsStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.RnaMods.txt");
            if (rnaModsStream != null)
            {
                using (var reader = new StreamReader(rnaModsStream))
                {
                    var rnaMods = ModificationLoader.ReadModsFromFile(reader,
                        new Dictionary<string, int>(), out var filteredWarnings);
                    AddModsToDictionary(allMods, rnaMods, "RNA");
                }
            }
        }
        catch (Exception ex)
        {
            System.Diagnostics.Debug.WriteLine($"Error loading RNA modifications: {ex.Message}");
            // Don't throw - RNA mods are optional
        }

        return allMods;
    }

    /// <summary>
    /// Helper method to add modifications to dictionary, handling duplicates
    /// </summary>
    private static void AddModsToDictionary(Dictionary<string, Modification> dict,
        IEnumerable<Modification> mods, string source)
    {
        foreach (var mod in mods)
        {
            if (!dict.ContainsKey(mod.IdWithMotif))
            {
                dict[mod.IdWithMotif] = mod;
            }
            else
            {
                // Log duplicate but prefer the first one loaded
                System.Diagnostics.Debug.WriteLine(
                    $"Duplicate modification from {source}: {mod.IdWithMotif}");
            }
        }
    }



    #endregion

    #region Public Methods

    /// <summary>
    /// Gets a modification by IdWithMotif or OriginalId
    /// </summary>
    public static Modification GetModification(string id, bool proteinOnly = false, bool rnaOnly = false)
    {
        if (proteinOnly)
        {
            if (AllKnownProteinMods.TryGetValue(id, out var mod))
                return mod;
            return AllProteinModsList.FirstOrDefault(m => m.OriginalId == id);
        }

        if (rnaOnly)
        {
            if (AllKnownRnaMods.TryGetValue(id, out var mod))
                return mod;
            return AllRnaModsList.FirstOrDefault(m => m.OriginalId == id);
        }

        // Search both
        if (AllModsKnown.TryGetValue(id, out var foundMod))
            return foundMod;

        return AllKnownMods.FirstOrDefault(m => m.OriginalId == id);
    }

    /// <summary>
    /// Gets all modifications matching a predicate
    /// </summary>
    public static IEnumerable<Modification> GetModifications(Func<Modification, bool> predicate,
        bool proteinOnly = false, bool rnaOnly = false)
    {
        if (proteinOnly)
            return AllProteinModsList.Where(predicate);

        if (rnaOnly)
            return AllRnaModsList.Where(predicate);

        return AllKnownMods.Where(predicate);
    }

    /// <summary>
    /// Adds or updates a custom modification (runtime only, not persisted)
    /// </summary>
    public static void AddOrUpdateModification(Modification modification, bool isRnaMod = false)
    {
        lock (_cacheLock)
        {
            if (isRnaMod)
            {
                AllKnownRnaMods[modification.IdWithMotif] = modification;
            }
            else
            {
                AllKnownProteinMods[modification.IdWithMotif] = modification;
            }
        }
    }

    /// <summary>
    /// Reloads all modifications from embedded resources
    /// Warning: This will clear any runtime-added modifications
    /// </summary>
    public static void ReloadModifications()
    {
        lock (_cacheLock)
        {
            // Force lazy re-initialization
            var lazyType = typeof(Lazy<Dictionary<string, Modification>>);
            var valueField = lazyType.GetField("m_value", BindingFlags.Instance | BindingFlags.NonPublic);

            valueField?.SetValue(_allKnownProteinMods, null);
            valueField?.SetValue(_allKnownRnaMods, null);

            var listLazyType = typeof(Lazy<List<Modification>>);
            var listValueField = listLazyType.GetField("m_value", BindingFlags.Instance | BindingFlags.NonPublic);

            listValueField?.SetValue(_allProteinModsList, null);
            listValueField?.SetValue(_allRnaModsList, null);
        }
    }

    #endregion
}
