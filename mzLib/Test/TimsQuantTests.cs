using FlashLFQ;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test
{ 
    [TestFixture]
    internal class TimsQuantTests
    {
        [Test]
        public static void TestFlashLfqTimsTof()
        {
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles");

            string dataFilePath = Path.Combine(testDataDirectory, "timsTOF_snippet.d");
            string psmFilePath = Path.Combine(testDataDirectory, "timsTOF_snippetSearch", "AllPSMs.psmtsv");
            string outputDirectory = Path.Combine(testDataDirectory, "timsTOF_snippetSearch", "flashTest");
            Directory.CreateDirectory(outputDirectory);

            SpectraFileInfo f1r1 = new SpectraFileInfo(dataFilePath, "one", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFilePath))
            {
                var split = line.Split(new char[] { '\t' });

                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = f1r1;

                string baseSequence = split[13];
                if (baseSequence.Contains("|")) continue;
                string fullSequence = split[14];
                double monoMass = double.Parse(split[23]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteinSubset = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new();

                if (fullSequence.Contains("[")) continue;

                foreach (var protein in proteinSubset)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: true,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: true,
                maxThreads: -1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), Path.Combine(outputDirectory, "bayesian.tsv"), true);

            var peaks = results.Peaks.Values.ToList();
            var peptides = results.PeptideModifiedSequences.Values.ToList();
            var proteins = results.ProteinGroups.Values.ToList();

            Assert.AreEqual(4, peaks[0].Count(m => m.IsMbrPeak == false));
            Assert.AreEqual(5, peaks[1].Count(m => m.IsMbrPeak == false));

            CollectionAssert.AreEquivalent(new string[] { "Q7KZF4", "Q7KZF4", "P52298", "Q15149", "Q15149" }, peaks[0].SelectMany(i => i.Identifications).Select(g => g.ProteinGroups.First()).Select(m => m.ProteinGroupName).ToArray());
            CollectionAssert.AreEquivalent(new string[] { "Q7KZF4", "P52298", "Q15149", "Q15149", "Q7KZF4", "Q7KZF4", "P52298" }, peaks[1].SelectMany(i => i.Identifications).Select(g => g.ProteinGroups.First()).Select(m => m.ProteinGroupName).ToArray());

            Assert.AreEqual(6, peptides.Count);
            CollectionAssert.AreEquivalent(new string[] { "Q7KZF4", "P52298", "Q15149", "Q15149", "Q7KZF4", "P52298" }, peptides.Select(g => g.ProteinGroups.First()).Select(m => m.ProteinGroupName).ToArray());

            Assert.AreEqual(3, proteins.Count);

            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).ToList();
            List<string> proteinsList = File.ReadAllLines(Path.Combine(outputDirectory, "proteins.tsv")).ToList();

            //check that all rows including header have the same number of elements
            Assert.AreEqual(1, peaksList.Select(l => l.Split('\t').Length).Distinct().ToList().Count);
            Assert.AreEqual(1, peptidesList.Select(l => l.Split('\t').Length).Distinct().ToList().Count);
            Assert.AreEqual(1, proteinsList.Select(l => l.Split('\t').Length).Distinct().ToList().Count);

            CollectionAssert.AreEquivalent(new string[] { "P52298", "Q15149", "Q7KZF4" }, proteins.Select(p => p.ProteinGroupName.ToArray()));

            //Directory.Delete(outputDirectory, true);
        }

    }
}
