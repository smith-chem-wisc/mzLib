using Omics.Modifications.IO;
using System.Reflection;

namespace Omics.Modifications;

public enum ModificationNamingConvention
{
    MetaMorpheus, 
    MetaMorpheus_Rna, 
    UniProt, 
    Unimod,
    Mixed
}

public static class Mods
{
    private static readonly Lazy<Dictionary<string, Modification>> _allKnownProteinMods;
    private static readonly Lazy<List<Modification>> _allProteinModsList;

    private static readonly Lazy<Dictionary<string, Modification>> _allKnownRnaMods;
    private static readonly Lazy<List<Modification>> _allRnaModsList;

    private static readonly object _cacheLock = new object();

    public static Dictionary<ModificationNamingConvention, List<Modification>> ModsByConvention { get; private set; }

    static Mods()
    {
        _allKnownProteinMods = new Lazy<Dictionary<string, Modification>>(LoadAllProteinModifications);
        _allProteinModsList = new Lazy<List<Modification>>(() => _allKnownProteinMods.Value.Values.ToList());
        _allKnownRnaMods = new Lazy<Dictionary<string, Modification>>(LoadAllRnaModifications);
        _allRnaModsList = new Lazy<List<Modification>>(() => _allKnownRnaMods.Value.Values.ToList());

        MetaMorpheusModifications = new List<Modification>();
        UnimodModifications = new List<Modification>();
        UniprotModifications = new List<Modification>();

        ModsByConvention = new Dictionary<ModificationNamingConvention, List<Modification>>
        {
            { ModificationNamingConvention.MetaMorpheus, MetaMorpheusModifications},
            { ModificationNamingConvention.MetaMorpheus_Rna, AllRnaModsList },
            { ModificationNamingConvention.UniProt, UniprotModifications},
            { ModificationNamingConvention.Unimod, UnimodModifications },
            { ModificationNamingConvention.Mixed, AllKnownMods }
        };
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

    public static List<Modification> UniprotModifications { get; private set; }
    public static List<Modification> MetaMorpheusModifications { get; private set; }
    public static List<Modification> UnimodModifications { get; private set; }

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
    
            UnimodModifications.Clear();
            var unimodMods = ModificationLoader.ReadModsFromUnimod(unimodStream!).ToList();

            UnimodModifications.AddRange(unimodMods);
            AddModsToDictionary(allMods, unimodMods, "Unimod");
      

            // 2. Load PSI-MOD and get formal charges
            var psiModStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.PSI-MOD.obo.xml");

            var psiModObo = ModificationLoader.LoadPsiMod(psiModStream!);
            Dictionary<string, int>  formalChargeDict = ModificationLoader.GetFormalChargesDictionary(psiModObo);
   

            // 3. Load UniProt ptmlist
            var uniprotStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.ptmlist.txt");

            UniprotModifications.Clear();
            using var uniProtReader = new StreamReader(uniprotStream!);
            var uniprotMods = ModificationLoader.ReadModsFromFile(uniProtReader, formalChargeDict,
                out _).ToList();

            UniprotModifications.AddRange(uniprotMods);
            AddModsToDictionary(allMods, uniprotMods, "UniProt");


            // 4. Load custom Mods.txt
            var modsTextStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.Mods.txt");

            MetaMorpheusModifications.Clear();
            using var modTextReader = new StreamReader(modsTextStream!);
            var modsTextMods = ModificationLoader.ReadModsFromFile(modTextReader, formalChargeDict,
                out _).ToList();

            MetaMorpheusModifications.AddRange(modsTextMods);
            AddModsToDictionary(allMods, modsTextMods, "MetaMorpheus ModsText");

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

            using var rnaModsReader = new StreamReader(rnaModsStream!);
            var rnaMods = ModificationLoader.ReadModsFromFile(rnaModsReader,
                new Dictionary<string, int>(), out _);
            AddModsToDictionary(allMods, rnaMods, "RNA");

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
    public static Modification? GetModification(string id, bool searchProteinMods = true, bool searchRnaMods = true)
    {
        if (!searchProteinMods && !searchRnaMods)
            throw new ArgumentException("At least one of searchProteinMods or searchRnaMods must be true.");

        if (searchProteinMods && !searchRnaMods)
        {
            if (AllKnownProteinMods.TryGetValue(id, out var mod))
                return mod;
            return AllProteinModsList.FirstOrDefault(m => m.IdWithMotif == id || m.OriginalId == id);
        }

        if (!searchProteinMods && searchRnaMods)
        {
            if (AllKnownRnaMods.TryGetValue(id, out var mod))
                return mod;
            return AllRnaModsList.FirstOrDefault(m => m.IdWithMotif == id || m.OriginalId == id);
        }

        // Search both
        if (AllModsKnown.TryGetValue(id, out var foundMod))
            return foundMod;

        return AllKnownMods.FirstOrDefault(m => m.IdWithMotif == id || m.OriginalId == id);
    }

    public static Modification? GetModification(string id, ModificationNamingConvention convention)
    {
        if (ModsByConvention.TryGetValue(convention, out var mods))
        {
            return mods.FirstOrDefault(m => m.IdWithMotif == id || m.OriginalId == id);
        }
        return null;
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

    #endregion
}
