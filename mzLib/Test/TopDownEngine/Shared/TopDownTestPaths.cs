using System.IO;
using NUnit.Framework;

namespace Test.TopDownEngine.Shared;

internal static class TopDownTestPaths
{
    internal static string Root => Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownEngine");
    internal static string Parsers => Path.Combine(Root, "Parsers");
    internal static string MetaMorpheus => Path.Combine(Parsers, "MetaMorpheus");
    internal static string ProSightPD => Path.Combine(Parsers, "ProSightPD");
    internal static string Consensus => Path.Combine(Root, "Consensus");
    internal static string Characterization => Path.Combine(Root, "Characterization");
    internal static string Simulation => Path.Combine(Root, "Simulation");
    internal static string Harness => Path.Combine(Root, "Harness");
    internal static string Regression => Path.Combine(Root, "Regression");

    internal static string Data => Path.GetFullPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "..", "..", "..", "..", "Data"));
    internal static string DataFiles => Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles");

    internal static string? FindExistingFile(params string[] relativePaths)
    {
        foreach (var candidate in relativePaths)
        {
            var dataPath = Path.Combine(Data, candidate);
            if (File.Exists(dataPath))
                return dataPath;

            var dataFilesPath = Path.Combine(DataFiles, candidate);
            if (File.Exists(dataFilesPath))
                return dataFilesPath;
        }

        return null;
    }

    internal static string GetRequiredFile(params string[] relativePaths) =>
        FindExistingFile(relativePaths) ?? throw new FileNotFoundException($"Could not find any of: {string.Join(", ", relativePaths)}");

    internal static void EnsureDirectory(string path)
    {
        if (!Directory.Exists(path))
            Directory.CreateDirectory(path);
    }
}
