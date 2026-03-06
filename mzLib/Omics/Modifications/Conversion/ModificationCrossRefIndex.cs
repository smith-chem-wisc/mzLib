using System;
using System.Collections.Concurrent;
using System.Collections.Immutable;
using System.Threading;

namespace Omics.Modifications.Conversion;

/// <summary>
/// Indexes modifications by database cross-reference identifiers to support
/// deterministic mapping between naming conventions.
/// </summary>
public sealed class ModificationCrossRefIndex
{
    private readonly IReadOnlyDictionary<string, ImmutableArray<Modification>> _databaseKeyIndex;
    private readonly IReadOnlyDictionary<string, ImmutableArray<Modification>> _accessionIndex;

    private static readonly object GlobalSync = new();
    private static ModificationCrossRefIndex? _global;
    private static int _globalVersion;

    private ModificationCrossRefIndex(
        IReadOnlyDictionary<string, ImmutableArray<Modification>> databaseKeyIndex,
        IReadOnlyDictionary<string, ImmutableArray<Modification>> accessionIndex)
    {
        _databaseKeyIndex = databaseKeyIndex;
        _accessionIndex = accessionIndex;
    }

    /// <summary>
    /// Gets the global index instance. Rebuilds automatically when the Mods registry version changes.
    /// </summary>
    public static ModificationCrossRefIndex Global
    {
        get
        {
            var cached = Volatile.Read(ref _global);
            if (cached != null && _globalVersion == Mods.RegistryVersion)
            {
                return cached;
            }

            lock (GlobalSync)
            {
                if (_global == null || _globalVersion != Mods.RegistryVersion)
                {
                    _global = Build(Mods.AllKnownMods);
                    _globalVersion = Mods.RegistryVersion;
                }

                return _global;
            }
        }
    }

    /// <summary>
    /// Forces a rebuild of the global index using the current Mods registry.
    /// </summary>
    public static void RefreshGlobal()
    {
        lock (GlobalSync)
        {
            _global = Build(Mods.AllKnownMods);
            _globalVersion = Mods.RegistryVersion;
        }
    }

    /// <summary>
    /// Gets all candidate modifications that share at least one cross-reference with the provided source modification.
    /// </summary>
    internal IReadOnlyList<Modification> GetCandidates(Modification source, Predicate<Modification> conventionFilter)
    {
        var candidates = new HashSet<Modification>();

        foreach (var key in EnumerateSourceKeys(source))
        {
            if (_databaseKeyIndex.TryGetValue(key, out var mods))
            {
                foreach (var mod in mods)
                {
                    if (conventionFilter(mod))
                    {
                        candidates.Add(mod);
                    }
                }
            }
        }

        if (!string.IsNullOrWhiteSpace(source.Accession))
        {
            var accessionKey = NormalizeValue(source.Accession);
            if (_accessionIndex.TryGetValue(accessionKey, out var accessionMods))
            {
                foreach (var mod in accessionMods)
                {
                    if (conventionFilter(mod))
                    {
                        candidates.Add(mod);
                    }
                }
            }
        }

        return candidates.ToList();
    }

    /// <summary>
    /// Gets all modifications registered under the specified database reference.
    /// </summary>
    public IReadOnlyList<Modification> GetByDatabaseId(string databaseName, string accession)
    {
        var key = NormalizeKey(databaseName, accession);
        if (_databaseKeyIndex.TryGetValue(key, out var mods))
        {
            return mods;
        }

        return Array.Empty<Modification>();
    }

    /// <summary>
    /// Gets all modifications registered with the provided accession identifier (e.g., PTM-0100).
    /// </summary>
    public IReadOnlyList<Modification> GetByAccession(string accession)
    {
        var normalized = NormalizeValue(accession);
        if (_accessionIndex.TryGetValue(normalized, out var mods))
        {
            return mods;
        }

        return Array.Empty<Modification>();
    }

    private static IEnumerable<string> EnumerateSourceKeys(Modification source)
    {
        if (source.DatabaseReference == null)
        {
            yield break;
        }

        foreach (var kvp in source.DatabaseReference)
        {
            if (string.IsNullOrWhiteSpace(kvp.Key) || kvp.Value == null)
            {
                continue;
            }

            var dbName = kvp.Key.Trim();
            foreach (var value in kvp.Value)
            {
                if (string.IsNullOrWhiteSpace(value))
                {
                    continue;
                }

                yield return NormalizeKey(dbName, value);
            }
        }
    }

    private static ModificationCrossRefIndex Build(IEnumerable<Modification> modifications)
    {
        var dbDict = new Dictionary<string, ImmutableArray<Modification>.Builder>(StringComparer.OrdinalIgnoreCase);
        var accessionDict = new Dictionary<string, ImmutableArray<Modification>.Builder>(StringComparer.OrdinalIgnoreCase);

        foreach (var mod in modifications)
        {
            if (mod.DatabaseReference != null)
            {
                foreach (var kvp in mod.DatabaseReference)
                {
                    if (string.IsNullOrWhiteSpace(kvp.Key) || kvp.Value == null)
                    {
                        continue;
                    }

                    foreach (var value in kvp.Value)
                    {
                        if (string.IsNullOrWhiteSpace(value))
                        {
                            continue;
                        }

                        var normalized = NormalizeKey(kvp.Key, value);
                        if (!dbDict.TryGetValue(normalized, out var builder))
                        {
                            builder = ImmutableArray.CreateBuilder<Modification>();
                            dbDict[normalized] = builder;
                        }

                        if (!builder.Contains(mod))
                        {
                            builder.Add(mod);
                        }
                    }
                }
            }

            if (!string.IsNullOrWhiteSpace(mod.Accession))
            {
                var normalizedAccession = NormalizeValue(mod.Accession);
                if (!accessionDict.TryGetValue(normalizedAccession, out var builder))
                {
                    builder = ImmutableArray.CreateBuilder<Modification>();
                    accessionDict[normalizedAccession] = builder;
                }

                if (!builder.Contains(mod))
                {
                    builder.Add(mod);
                }
            }
        }

        var dbIndex = dbDict.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.ToImmutable(), StringComparer.OrdinalIgnoreCase);
        var accessionIndex = accessionDict.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.ToImmutable(), StringComparer.OrdinalIgnoreCase);

        return new ModificationCrossRefIndex(dbIndex, accessionIndex);
    }

    private static string NormalizeKey(string databaseName, string accession)
    {
        return $"{NormalizeValue(databaseName)}:{NormalizeValue(accession)}";
    }

    private static string NormalizeValue(string value)
    {
        return value.Trim().ToUpperInvariant();
    }
}
