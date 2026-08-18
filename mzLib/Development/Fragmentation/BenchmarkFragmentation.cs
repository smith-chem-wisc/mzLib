// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (TestFragments.cs) is part of Proteomics.
//
// Proteomics is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Proteomics is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Proteomics. If not, see <http://www.gnu.org/licenses/>.

using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Development
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class BenchmarkFragmentation
    {
        
        [Test]
        public static void Benchmark_ParallelFragmentation()
        {
            // Load proteins from the cRAP database with reverse decoys to maximize peptide count
            //var dbPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cRAP_databaseGPTMD.xml");

            // This path points to a human database where I have run GPTMD with a large number of variable modifications, which should yield a very large number of peptides and fragments to maximize the fragmentation time and thus the potential speedup from parallelization. Adjust the path as needed to point to a suitable test database on your machine.
            // I'm specifically interested in benchmarking performance with neutral losses, and a standard uniprot xml does not contain neutral losses
            var dbPath = @"D:\Kelly_ALS_motor_nueron_dataset\MM1p1p4_GPTMD_Search_Carboxymethyl_Carbamido_2\Task1-GPTMDTask\uniprotkb_Human_AND_model_organism_9606_2025_03_19GPTMD.xml";

            var loadSw = Stopwatch.StartNew();
            var proteins = ProteinDbLoader.LoadProteinXML(dbPath, true, DecoyType.Reverse, Mods.AllKnownMods, false, null, out _, maxHeterozygousVariants: 0);
            loadSw.Stop();
            var loadElapsed = loadSw.Elapsed;

            // Digest all proteins into peptides
            var digestionParams = new DigestionParams();
            var peptides = proteins
                .SelectMany(p => p.Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                .ToList();

            // Warm up: trigger JIT compilation and lazy initialization before timing
            var warmupProducts = new List<Product>();
            peptides[0].Fragment(DissociationType.HCD, FragmentationTerminus.Both, warmupProducts);

            // Serial benchmark
            long serialFragmentCount = 0;
            var serialProducts = new List<Product>();
            var sw = Stopwatch.StartNew();
            foreach (var peptide in peptides)
            {
                peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, serialProducts);
                serialFragmentCount += serialProducts.Count;
            }
            sw.Stop();
            var serialElapsed = sw.Elapsed;

            // Parallel benchmark: each thread owns its List<Product> to avoid contention
            long parallelFragmentCount = 0;
            sw.Restart();
            Parallel.ForEach(
                peptides,
                () => new List<Product>(),
                (peptide, _, localProducts) =>
                {
                    peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, localProducts);
                    Interlocked.Add(ref parallelFragmentCount, localProducts.Count);
                    return localProducts;
                },
                _ => { });
            sw.Stop();
            var parallelElapsed = sw.Elapsed;

            TestContext.Out.WriteLine($"Proteins loaded:          {proteins.Count}");
            TestContext.Out.WriteLine($"Loading elapsed:          {loadElapsed.TotalSeconds:F3} s");
            TestContext.Out.WriteLine($"Peptides digested:        {peptides.Count}");
            TestContext.Out.WriteLine($"Serial fragment count:    {serialFragmentCount}");
            TestContext.Out.WriteLine($"Serial elapsed:           {serialElapsed.TotalSeconds:F3} s");
            TestContext.Out.WriteLine($"Parallel fragment count:  {parallelFragmentCount}");
            TestContext.Out.WriteLine($"Parallel elapsed:         {parallelElapsed.TotalSeconds:F3} s");
            TestContext.Out.WriteLine($"Speedup:                  {serialElapsed.TotalSeconds / parallelElapsed.TotalSeconds:F2}x");

            Assert.That(serialFragmentCount, Is.GreaterThan(0));
            Assert.That(serialFragmentCount, Is.EqualTo(parallelFragmentCount));
        }
    }
}