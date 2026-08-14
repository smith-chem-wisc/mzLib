using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection.Metadata;
using System.Reflection.PortableExecutable;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ReadersInfrastructure
{
    /// <summary>
    /// A managed assembly stamped for one architecture cannot be loaded by CoreCLR on another, so using a
    /// type from one makes the loader throw BadImageFormatException far from the call site. Referencing it
    /// is harmless on its own: the C# compiler omits references that no code uses, and the CLR resolves
    /// assemblies lazily, so an unused reference never loads.
    ///
    /// That is the situation Thermo's MassPrecisionEstimator was in - referenced by Readers.csproj, stamped
    /// AMD64, used by nothing, and one call away from breaking .raw reading on arm64. The reference is gone;
    /// these tests check the property rather than the name, so a future architecture-neutral build of any
    /// such assembly is free to be used.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestNoArchitectureSpecificAssemblyReferences
    {
        /// <summary>
        /// An architecture-neutral managed assembly is IL-only and leaves the PE machine field as I386,
        /// which the runtime reads as "any architecture". Anything else - AMD64 or ARM64 stamping, or a
        /// mixed-mode assembly carrying native code - loads only where it was built for.
        /// </summary>
        private static bool IsArchitectureNeutral(PEReader peReader) =>
            peReader.PEHeaders.CorHeader.Flags.HasFlag(CorFlags.ILOnly)
            && peReader.PEHeaders.CoffHeader.Machine == Machine.I386;

        private static string Stamp(PEReader peReader) =>
            $"machine={peReader.PEHeaders.CoffHeader.Machine}, " +
            $"ilOnly={peReader.PEHeaders.CorHeader.Flags.HasFlag(CorFlags.ILOnly)}";

        /// <summary>
        /// Every file in the build output that has managed metadata, keyed by assembly simple name. Native
        /// libraries (isodeclib, timsdata, baf2sql) have no metadata and are skipped: they are loaded by
        /// P/Invoke, which fails by a different route this test cannot see.
        /// </summary>
        private static Dictionary<string, string> ManagedAssembliesByName()
        {
            Dictionary<string, string> byName = new(StringComparer.OrdinalIgnoreCase);

            foreach (string path in Directory
                .EnumerateFiles(TestContext.CurrentContext.TestDirectory, "*.dll", SearchOption.AllDirectories)
                .OrderBy(path => path))
            {
                try
                {
                    using FileStream stream = File.OpenRead(path);
                    using PEReader peReader = new(stream);
                    if (!peReader.HasMetadata)
                        continue;
                }
                catch (BadImageFormatException)
                {
                    continue;
                }

                byName.TryAdd(Path.GetFileNameWithoutExtension(path), path);
            }

            return byName;
        }

        private static List<string> ReferencedAssemblyNames(string path)
        {
            using FileStream stream = File.OpenRead(path);
            using PEReader peReader = new(stream);
            MetadataReader metadata = peReader.GetMetadataReader();
            return metadata.AssemblyReferences
                .Select(handle => metadata.GetString(metadata.GetAssemblyReference(handle).Name))
                .ToList();
        }

        [Test]
        public void NoAssemblyInTheBuildOutputReferencesAnArchitectureSpecificAssembly()
        {
            Dictionary<string, string> byName = ManagedAssembliesByName();
            List<string> offenders = new();

            foreach ((string name, string path) in byName)
            {
                foreach (string referenceName in ReferencedAssemblyNames(path))
                {
                    // a reference we cannot find in the output cannot be classified; it would fail to load
                    // for its own reasons long before architecture mattered
                    if (!byName.TryGetValue(referenceName, out string referencePath))
                        continue;

                    using FileStream stream = File.OpenRead(referencePath);
                    using PEReader referenced = new(stream);
                    if (!IsArchitectureNeutral(referenced))
                        offenders.Add($"{name}.dll -> {referenceName} ({Stamp(referenced)})");
                }
            }

            // guards against passing because the scan found nothing to examine
            Assert.That(byName, Has.Count.GreaterThan(10),
                "expected to scan the build output; the directory looks wrong or empty");

            Assert.That(offenders, Is.Empty,
                "these references are to assemblies stamped for a single architecture, so using them breaks " +
                "on every other architecture:\n  " + string.Join("\n  ", offenders));
        }

        /// <summary>
        /// Positive control for the classifier. Without a genuinely architecture-stamped input, the scan
        /// above would pass whether or not its logic worked. ThermoFisher.CommonCore.MassPrecisionEstimator
        /// is kept in the test output purely as that input - nothing references it, and it is no longer
        /// shipped in the package.
        /// </summary>
        [Test]
        public void TheClassifierRecognisesAnArchitectureStampedAssembly()
        {
            string x64Assembly = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "FileReadingTests", "ArchitectureGuard", "ThermoFisher.CommonCore.MassPrecisionEstimator.dll");

            Assert.That(File.Exists(x64Assembly), Is.True, $"missing test input: {x64Assembly}");

            using FileStream stream = File.OpenRead(x64Assembly);
            using PEReader peReader = new(stream);

            Assert.That(peReader.HasMetadata, Is.True, "expected a managed assembly");
            Assert.That(IsArchitectureNeutral(peReader), Is.False,
                $"this assembly is stamped AMD64 and must be classified as architecture-specific ({Stamp(peReader)})");
        }

        /// <summary>
        /// Positive control for reference resolution: a reference that is genuinely in use has to be found
        /// in the output and classified, otherwise the scan would silently skip everything.
        /// </summary>
        [Test]
        public void TheScanResolvesAndClassifiesAReferenceThatIsInUse()
        {
            string readersPath = typeof(ThermoRawFileReader).Assembly.Location;
            List<string> references = ReferencedAssemblyNames(readersPath);

            Assert.That(references, Does.Contain("ThermoFisher.CommonCore.Data"),
                "Readers uses types from ThermoFisher.CommonCore.Data, so that reference must be visible");

            Dictionary<string, string> byName = ManagedAssembliesByName();
            Assert.That(byName.ContainsKey("ThermoFisher.CommonCore.Data"), Is.True,
                "the reference must resolve to a file in the output for the scan to classify it");

            using FileStream stream = File.OpenRead(byName["ThermoFisher.CommonCore.Data"]);
            using PEReader peReader = new(stream);
            Assert.That(IsArchitectureNeutral(peReader), Is.True,
                "ThermoFisher.CommonCore.Data is architecture-neutral and must not be flagged");
        }
    }
}
