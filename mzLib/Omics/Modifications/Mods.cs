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
    private static readonly object _cacheLock = new object();

    public static Dictionary<ModificationNamingConvention, List<Modification>> ModsByConvention { get; private set; }

    static Mods()
    {
        LoadAllProteinModifications();
        AllProteinModsList = UnimodModifications.Concat(UniprotModifications).Concat(MetaMorpheusProteinModifications).ToList();
        AllKnownProteinModsDictionary = AllProteinModsList
            .DistinctBy(m => m.IdWithMotif)
            .ToDictionary(m => m.IdWithMotif);

        LoadAllRnaModifications();
        AllRnaModsList = MetaMorpheusRnaModifications.ToList();
        AllKnownRnaModsDictionary = AllRnaModsList
            .DistinctBy(m => m.IdWithMotif)
            .ToDictionary(m => m.IdWithMotif);

        // Combine protein and RNA mods, with Protein mods taking precedence in case of conflicts
        AllKnownMods = AllProteinModsList.Concat(AllRnaModsList).ToList();
        AllModsKnownDictionary = new Dictionary<string, Modification>(AllKnownRnaModsDictionary);
        foreach (var kvp in AllKnownProteinModsDictionary)
        {
            AllModsKnownDictionary[kvp.Key] = kvp.Value;
        }

        ModsByConvention = new Dictionary<ModificationNamingConvention, List<Modification>>
        {
            { ModificationNamingConvention.MetaMorpheus, MetaMorpheusProteinModifications},
            { ModificationNamingConvention.MetaMorpheus_Rna, MetaMorpheusRnaModifications },
            { ModificationNamingConvention.UniProt, UniprotModifications},
            { ModificationNamingConvention.Unimod, UnimodModifications },
            { ModificationNamingConvention.Mixed, AllKnownMods }
        };
    }

    #region Public Properties
    public static List<Modification> UniprotModifications { get; private set; } = [];
    public static List<Modification> MetaMorpheusProteinModifications { get; private set; } = [];
    public static List<Modification> UnimodModifications { get; private set; } = [];

    /// <summary>
    /// All known protein modifications indexed by IdWithMotif
    /// </summary>
    public static Dictionary<string, Modification> AllKnownProteinModsDictionary { get; }

    /// <summary>
    /// All known protein modifications as a list
    /// </summary>
    public static List<Modification> AllProteinModsList { get; }



    public static List<Modification> MetaMorpheusRnaModifications { get; private set; } = [];

    /// <summary>
    /// All known RNA modifications indexed by IdWithMotif
    /// </summary>
    public static Dictionary<string, Modification> AllKnownRnaModsDictionary { get; }

    /// <summary>
    /// All known RNA modifications as a list
    /// </summary>
    public static List<Modification> AllRnaModsList { get;  }

    /// <summary>
    /// Combined list of all known modifications (protein + RNA)
    /// </summary>
    public static List<Modification> AllKnownMods { get; }

    /// <summary>
    /// Combined dictionary of all known modifications (protein + RNA)
    /// Protein mods take precedence in case of conflicts
    /// </summary>
    public static Dictionary<string, Modification> AllModsKnownDictionary { get; }


    #endregion

    #region Loading Methods

    /// <summary>
    /// Loads Protein modifications from embedded resources and custom sources
    /// </summary>
    private static void LoadAllProteinModifications()
    {
        var assembly = Assembly.GetExecutingAssembly();
        var assemblyName = assembly.GetName().Name;

        // 1. Load Unimod
        var unimodStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.unimod.xml");
    
        UnimodModifications = ModificationLoader.ReadModsFromUnimod(unimodStream!).ToList();   

        // 2. Load PSI-MOD and get formal charges
        var psiModStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.PSI-MOD.obo.xml");

        var psiModObo = ModificationLoader.LoadPsiMod(psiModStream!);
        Dictionary<string, int>  formalChargeDict = ModificationLoader.GetFormalChargesDictionary(psiModObo);
   

        // 3. Load UniProt ptmlist
        var uniprotStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.ptmlist.txt");
        using var uniProtReader = new StreamReader(uniprotStream!);
        UniprotModifications = ModificationLoader.ReadModsFromFile(uniProtReader, formalChargeDict,
            out _).ToList();

        // 4. Load custom Mods.txt
        var modsTextStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.Mods.txt");
        using var modTextReader = new StreamReader(modsTextStream!);
        MetaMorpheusProteinModifications = ModificationLoader.ReadModsFromFile(modTextReader, formalChargeDict,
            out _).ToList();
    }

    /// <summary>
    /// Loads RNA modifications from embedded resources and custom sources
    /// </summary>
    private static void LoadAllRnaModifications()
    {
        var assembly = Assembly.GetExecutingAssembly();
        var assemblyName = assembly.GetName().Name;
        var rnaModsStream = assembly.GetManifestResourceStream($"{assemblyName}.Resources.RnaMods.txt");

        using var rnaModsReader = new StreamReader(rnaModsStream!);
        MetaMorpheusRnaModifications = ModificationLoader.ReadModsFromFile(rnaModsReader, new Dictionary<string, int>(), out _).ToList();
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
            if (AllKnownProteinModsDictionary.TryGetValue(id, out var mod))
                return mod;
            return AllProteinModsList.FirstOrDefault(m => m.IdWithMotif == id);
        }

        if (!searchProteinMods && searchRnaMods)
        {
            if (AllKnownRnaModsDictionary.TryGetValue(id, out var mod))
                return mod;
            return AllRnaModsList.FirstOrDefault(m => m.IdWithMotif == id);
        }

        // Search both
        if (AllModsKnownDictionary.TryGetValue(id, out var foundMod))
            return foundMod;

        return AllKnownMods.FirstOrDefault(m => m.IdWithMotif == id);
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
                AllKnownRnaModsDictionary[modification.IdWithMotif] = modification;
            }
            else
            {
                AllKnownProteinModsDictionary[modification.IdWithMotif] = modification;
            }
        }
    }

    #endregion
}
