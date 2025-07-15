using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Digestion;
using Omics.Modifications;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestSeqCoverage
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void MultipleProteaseSelectionTest()
        {
            Protein ParentProtein = new Protein("MOAT", "accession1");
            var motifList = DigestionMotif.ParseDigestionMotifsFromString("O|,|T");
            var protease = new Protease("TestProtease1", CleavageSpecificity.Full, null, null, motifList);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams multiProtease = new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList = ParentProtein.Digest(multiProtease, new List<Modification>(), new List<Modification>()).ToList();
            var sequences = digestedList.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 3);
            Assert.That(sequences.Contains("MO"));
            Assert.That(sequences.Contains("A"));
            Assert.That(sequences.Contains("T"));
        }

        [Test]
        public static void MultipleProteaseSelectionTestMissedCleavage()
        {
            Protein ParentProtein = new Protein("MOAT", "accession1");

            var motifList = DigestionMotif.ParseDigestionMotifsFromString("O|,|T");           
            var protease = new Protease("TestProtease2", CleavageSpecificity.Full, null, null, motifList);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams multiProtease = new DigestionParams(protease: protease.Name, maxMissedCleavages: 1, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList = ParentProtein.Digest(multiProtease, new List<Modification>(), new List<Modification>()).ToList();
            var sequences = digestedList.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 5);
            Assert.That(sequences.Contains("MOA"));
            Assert.That(sequences.Contains("AT"));
            Assert.That(sequences.Contains("MO"));
            Assert.That(sequences.Contains("A"));
            Assert.That(sequences.Contains("T"));
        }

        [Test]
        public static void MultipleProteaseSelectionTestPreventCleavage()
        {
            Protein ParentProtein = new Protein("MOAT", "accession1");

            var motifList = DigestionMotif.ParseDigestionMotifsFromString("O[A]|,|T");
            var protease = new Protease("TestProtease3", CleavageSpecificity.Full, null, null, motifList);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams multiProtease = new DigestionParams(protease: protease.Name, maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList = ParentProtein.Digest(multiProtease, new List<Modification>(), new List<Modification>()).ToList();
            var sequences = digestedList.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 2);
            Assert.That(sequences.Contains("MOA"));
            Assert.That(sequences.Contains("T"));
        }

        [Test]
        public static void ReadCustomFile()
        {
            Protein ParentProtein = new Protein("OKAREDY", "accession1");
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DoubleProtease.tsv");
            Assert.That(File.Exists(path));

            var proteaseDict = ProteaseDictionary.LoadProteaseDictionary(path);
            Assert.That(proteaseDict.ContainsKey("Test1"));
            Assert.That(proteaseDict.ContainsKey("Test2"));
            Assert.That(proteaseDict.ContainsKey("Test3"));
            ProteaseDictionary.Dictionary.Add("Test1", proteaseDict["Test1"]);

            DigestionParams multiProtease1 = new DigestionParams(protease: "Test1", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList1 = ParentProtein.Digest(multiProtease1, new List<Modification>(), new List<Modification>()).ToList();
            ProteaseDictionary.Dictionary.Remove("Test1");

            var sequences = digestedList1.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences.Count == 3);
            Assert.That(sequences.Contains("OK"));
            Assert.That(sequences.Contains("A"));
            Assert.That(sequences.Contains("REDY"));

            ProteaseDictionary.Dictionary.Add("Test2", proteaseDict["Test2"]);
            DigestionParams multiProtease2 = new DigestionParams(protease: "Test2", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList2 = ParentProtein.Digest(multiProtease2, new List<Modification>(), new List<Modification>()).ToList();
            ProteaseDictionary.Dictionary.Remove("Test2");
            var sequences2 = digestedList2.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences2.Count == 3);
            Assert.That(sequences2.Contains("OK"));
            Assert.That(sequences2.Contains("ARED"));
            Assert.That(sequences2.Contains("Y"));

            ProteaseDictionary.Dictionary.Add("Test3", proteaseDict["Test3"]);
            DigestionParams multiProtease3 = new DigestionParams(protease: "Test3", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var digestedList3 = ParentProtein.Digest(multiProtease3, new List<Modification>(), new List<Modification>()).ToList();
            ProteaseDictionary.Dictionary.Remove("Test3");
            var sequences3 = digestedList3.Select(p => p.BaseSequence).ToList();
            Assert.That(sequences3.Count == 2);
            Assert.That(sequences3.Contains("OK"));
            Assert.That(sequences3.Contains("AREDY"));
        }
    }
}