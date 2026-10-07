using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// The modification data files are described in two places that have to agree and never meet:
    /// mzLib.nuspec decides where they land inside the .nupkg, and nuget/mzLib.targets decides where
    /// a consuming build looks for them. Nothing at pack time compares the two, so if either side is
    /// renamed the pack still succeeds and consumers who opt in silently get nothing.
    ///
    /// These read the two files as text rather than packing, so they run anywhere and fail fast.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class NuGetPackagingTests
    {
        /// <summary>
        /// Walks up from the test output directory to the folder holding mzLib.nuspec.
        /// </summary>
        private static string RepositoryRoot()
        {
            var dir = new DirectoryInfo(TestContext.CurrentContext.TestDirectory);
            while (dir != null && !File.Exists(Path.Combine(dir.FullName, "mzLib.nuspec")))
            {
                dir = dir.Parent;
            }

            Assert.That(dir, Is.Not.Null, "could not find the directory containing mzLib.nuspec");
            return dir.FullName;
        }

        private static string NuspecPath() => Path.Combine(RepositoryRoot(), "mzLib.nuspec");

        private static string TargetsPath() => Path.Combine(RepositoryRoot(), "nuget", "mzLib.targets");

        /// <summary>
        /// Every modification data file the nuspec ships must land under the root the targets file
        /// reads from, in a subfolder the targets file actually globs.
        /// </summary>
        [Test]
        public static void ModificationDataLayoutAgreesBetweenNuspecAndTargets()
        {
            var nuspec = XDocument.Load(NuspecPath());
            var targets = File.ReadAllText(TargetsPath());

            string[] shippedTargets = nuspec.Descendants()
                .Where(e => e.Name.LocalName == "file")
                .Select(e => (string)e.Attribute("target"))
                .Where(t => t != null && t.StartsWith("modificationData", StringComparison.Ordinal))
                .Distinct()
                .ToArray();

            Assert.That(shippedTargets, Is.Not.Empty,
                "the nuspec ships no modificationData files; if that is deliberate, this fixture should go too");

            // <MzLibModificationDataRoot>$(MSBuildThisFileDirectory)..\modificationData\</...>
            Match root = Regex.Match(targets, @"<MzLibModificationDataRoot>[^<]*?\\(?<folder>[^\\<]+)\\</MzLibModificationDataRoot>");
            Assert.That(root.Success, Is.True, "could not read MzLibModificationDataRoot out of mzLib.targets");
            Assert.That(root.Groups["folder"].Value, Is.EqualTo("modificationData"),
                "the targets file reads from a different folder than the nuspec writes to");

            // <None Include="$(MzLibModificationDataRoot)Mods\*" ... />
            var globbed = Regex.Matches(targets, @"\$\(MzLibModificationDataRoot\)(?<sub>[^\\""]+)\\\*")
                .Select(m => m.Groups["sub"].Value)
                .ToHashSet(StringComparer.Ordinal);

            Assert.That(globbed, Is.Not.Empty, "mzLib.targets globs no subfolders under the data root");

            var shippedSubfolders = shippedTargets
                .Select(t => t.Split('\\').Skip(1).FirstOrDefault())
                .Where(s => !string.IsNullOrEmpty(s))
                .Distinct(StringComparer.Ordinal)
                .ToArray();

            Assert.That(shippedSubfolders, Is.EquivalentTo(globbed),
                "the nuspec ships modificationData subfolders that mzLib.targets does not copy, or vice versa");
        }

        /// <summary>
        /// A nuspec file entry whose source is missing fails the pack, so a stale src is a broken
        /// release rather than a broken consumer. Cheaper to catch here.
        /// </summary>
        [Test]
        public static void EveryModificationDataFileTheNuspecShipsExists()
        {
            string root = RepositoryRoot();
            var missing = new List<string>();

            foreach (var file in XDocument.Load(NuspecPath()).Descendants()
                         .Where(e => e.Name.LocalName == "file"))
            {
                string target = (string)file.Attribute("target") ?? "";
                if (!target.StartsWith("modificationData", StringComparison.Ordinal)) continue;

                string src = ((string)file.Attribute("src") ?? "").Replace('\\', Path.DirectorySeparatorChar);
                if (!File.Exists(Path.Combine(root, src))) missing.Add(src);
            }

            Assert.That(missing, Is.Empty, "nuspec lists modification data files that do not exist");
        }
    }
}
