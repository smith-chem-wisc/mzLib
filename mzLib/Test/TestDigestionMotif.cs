using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestDigestionMotif
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
        public static void TestParseProtease()
        {
            var argn = DigestionMotif.ParseDigestionMotifsFromString("|D");
            Assert.AreEqual(argn.Count, 1);
            var c = argn[0];
            Assert.AreEqual(c.InducingCleavage, "D");
            Assert.AreEqual(c.PreventingCleavage, null);
            Assert.AreEqual(c.CutIndex, 0);
            var chymotrypsin = DigestionMotif.ParseDigestionMotifsFromString("F[P]|,W[P]|,Y[P]|");
            Assert.AreEqual(chymotrypsin.Count, 3);
        }
        [Test]
        public static void TestAllProteasesAgainstAllPeptides()
        {
            // File paths
            string fastaPath = @"E:\Projects\WisperBio\20251201_JamesImmuto\SevenProtein.fasta";
            string observedPeptidesPath = @"E:\Projects\WisperBio\20251201_JamesImmuto\2026-01-28-13-53-36\Task1-SearchTask\observedPeptides.tsv";
            string outputFolder = @"E:\Projects\WisperBio\20251201_JamesImmuto\2026-01-28-13-53-36\Task1-SearchTask";
            string poorlyObservedPath = Path.Combine(outputFolder, "poorlyObservedRegions.tsv");
            string peptidesNotSeenPath = Path.Combine(outputFolder, "peptidesNotSeen.tsv");

            // Load proteins from FASTA
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(fastaPath, true, DecoyType.None, false,
                out var dbErrors);
            var proteinDict = proteins.ToDictionary(p => p.Accession, p => p);

            // Digestion parameters for each protease
            string[] proteases = { "trypsin|P", "chymotrypsin|P", "elastase|P" };
            var emptyMods = new List<Modification>();

            // Digest proteins with all proteases and collect peptides with their positions
            var theoreticalPeptidesByProtein = new Dictionary<string, List<(int Start, int End, string Sequence, int MissedCleavages)>>();

            foreach (var protein in proteins)
            {
                theoreticalPeptidesByProtein[protein.Accession] = new List<(int, int, string, int)>();

                foreach (var proteaseName in proteases)
                {
                    var digestionParams = new DigestionParams(
                        protease: proteaseName,
                        maxMissedCleavages: 1,
                        minPeptideLength: 6,
                        searchModeType: CleavageSpecificity.Full);

                    var peptides = protein.Digest(digestionParams, emptyMods, emptyMods);

                    foreach (var peptide in peptides)
                    {
                        theoreticalPeptidesByProtein[protein.Accession].Add((
                            peptide.OneBasedStartResidueInProtein,
                            peptide.OneBasedEndResidueInProtein,
                            peptide.BaseSequence,
                            peptide.MissedCleavages));
                    }
                }

                // Remove duplicates (same peptide from different proteases)
                theoreticalPeptidesByProtein[protein.Accession] = theoreticalPeptidesByProtein[protein.Accession]
                    .Distinct()
                    .ToList();
            }

            // Read observed peptides (accession in col 1, peptide in col 2)
            var observedPeptidesByProtein = new Dictionary<string, List<(int Start, int End, string Sequence)>>();
            var observedPeptideSet = new HashSet<(string Accession, string Sequence)>();

            foreach (var line in File.ReadLines(observedPeptidesPath).Skip(1)) // Skip header
            {
                var parts = line.Split('\t');
                if (parts.Length < 2) continue;

                string accession = parts[0].Trim();
                string peptideSequence = parts[1].Trim();

                if (!proteinDict.TryGetValue(accession, out var protein)) continue;

                // Track observed peptides for comparison
                observedPeptideSet.Add((accession, peptideSequence));

                // Find peptide position in protein by string matching
                int index = protein.BaseSequence.IndexOf(peptideSequence, StringComparison.Ordinal);
                if (index < 0) continue;

                int oneBasedStart = index + 1;
                int oneBasedEnd = index + peptideSequence.Length;

                if (!observedPeptidesByProtein.ContainsKey(accession))
                    observedPeptidesByProtein[accession] = new List<(int, int, string)>();

                observedPeptidesByProtein[accession].Add((oneBasedStart, oneBasedEnd, peptideSequence));
            }

            // Write poorly observed regions
            using (var writer = new StreamWriter(poorlyObservedPath))
            {
                writer.WriteLine("Accession\tRegion\tSequence");

                foreach (var protein in proteins)
                {
                    string accession = protein.Accession;
                    int proteinLength = protein.BaseSequence.Length;

                    // Initialize coverage array (count of peptides covering each position)
                    int[] coverage = new int[proteinLength];

                    // Add coverage from observed peptides
                    if (observedPeptidesByProtein.TryGetValue(accession, out var observedPeptides))
                    {
                        foreach (var (start, end, _) in observedPeptides)
                        {
                            for (int i = start - 1; i < end; i++) // Convert to 0-based
                            {
                                coverage[i]++;
                            }
                        }
                    }

                    // Find contiguous regions with coverage <= 1
                    var lowCoverageRegions = new List<(int Start, int End)>();
                    int? regionStart = null;

                    for (int i = 0; i < proteinLength; i++)
                    {
                        if (coverage[i] <= 1)
                        {
                            if (regionStart == null)
                                regionStart = i + 1; // Convert to 1-based
                        }
                        else
                        {
                            if (regionStart != null)
                            {
                                lowCoverageRegions.Add((regionStart.Value, i)); // End is 1-based
                                regionStart = null;
                            }
                        }
                    }

                    // Close final region if it extends to the end
                    if (regionStart != null)
                    {
                        lowCoverageRegions.Add((regionStart.Value, proteinLength));
                    }

                    // Write results to file
                    foreach (var (start, end) in lowCoverageRegions)
                    {
                        string segment = protein.BaseSequence.Substring(start - 1, end - start + 1);
                        writer.WriteLine($"{accession}\t{start}-{end}\t{segment}");
                    }
                }
            }

            // Write peptides not seen (theoretical peptides not found in observed)
            using (var writer = new StreamWriter(peptidesNotSeenPath))
            {
                writer.WriteLine("Accession\tStart\tEnd\tSequence\tMissedCleavages");

                foreach (var protein in proteins)
                {
                    string accession = protein.Accession;

                    if (!theoreticalPeptidesByProtein.TryGetValue(accession, out var theoreticalPeptides))
                        continue;

                    foreach (var (start, end, sequence, missedCleavages) in theoreticalPeptides)
                    {
                        // Check if this peptide was observed
                        if (!observedPeptideSet.Contains((accession, sequence)))
                        {
                            writer.WriteLine($"{accession}\t{start}\t{end}\t{sequence}\t{missedCleavages}");
                        }
                    }
                }
            }

            // Console summary
            Console.WriteLine($"Poorly observed regions written to: {poorlyObservedPath}");
            Console.WriteLine($"Peptides not seen written to: {peptidesNotSeenPath}");
            Console.WriteLine("\n=== SUMMARY ===");

            int totalTheoretical = 0;
            int totalObserved = 0;
            int totalNotSeen = 0;

            foreach (var protein in proteins)
            {
                string accession = protein.Accession;
                var theoretical = theoreticalPeptidesByProtein.GetValueOrDefault(accession) ?? new List<(int, int, string, int)>();
                int theoreticalCount = theoretical.Count;
                int observedCount = observedPeptidesByProtein.GetValueOrDefault(accession)?.Count ?? 0;
                int notSeenCount = theoretical.Count(t => !observedPeptideSet.Contains((accession, t.Sequence)));

                totalTheoretical += theoreticalCount;
                totalObserved += observedCount;
                totalNotSeen += notSeenCount;

                // Calculate overall coverage percentage
                int[] coverage = new int[protein.BaseSequence.Length];
                if (observedPeptidesByProtein.TryGetValue(accession, out var observed))
                {
                    foreach (var (start, end, _) in observed)
                    {
                        for (int i = start - 1; i < end; i++)
                            coverage[i]++;
                    }
                }

                int coveredResidues = coverage.Count(c => c > 1);
                double coveragePercent = 100.0 * coveredResidues / protein.BaseSequence.Length;

                Console.WriteLine($"{accession}: {observedCount} observed / {theoreticalCount} theoretical peptides, " +
                                  $"{notSeenCount} not seen, {coveragePercent:F1}% residues with >1 coverage");
            }

            Console.WriteLine($"\nTOTAL: {totalObserved} observed / {totalTheoretical} theoretical, {totalNotSeen} peptides not seen");

            // Basic assertions
            Assert.IsTrue(proteins.Count > 0, "No proteins loaded from FASTA");
            Assert.IsTrue(observedPeptidesByProtein.Count > 0, "No observed peptides loaded");
            Assert.IsTrue(File.Exists(poorlyObservedPath), "Poorly observed regions file was not created");
            Assert.IsTrue(File.Exists(peptidesNotSeenPath), "Peptides not seen file was not created");
        }

        [Test]
        public static void TestBasicProtease1()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams(minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("PROTEIN", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "PR");
            Assert.AreEqual(last, "OTEIN");
        }

        [Test]
        public static void TestBasicProtease2()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("Lys-C|P", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("MKPKPKPMKA", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual("MKPKPKPMK", first);
            Assert.AreEqual("A", last);
        }

        [Test]
        public static void TestWildCardExclusion()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("RX{P}|");
            Protease multiletter = new Protease("multiletter", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(multiletter.Name, multiletter);
            DigestionParams myDigestionParams = new DigestionParams("multiletter", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("PROPRPPM", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "PRO");
            Assert.AreEqual(last, "PRPPM");
        }

        [Test]
        public static void TestMultiLetterNTerm()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("|AAA");
            Protease multiletter = new Protease("multi-custom", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(multiletter.Name, multiletter);
            DigestionParams myDigestionParams = new DigestionParams("multi-custom", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("FAAAMAAM", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "F");
            Assert.AreEqual(last, "AAAMAAM");
        }

        [Test]
        public static void TestMultiLetterProtease()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("collagenase", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("ABCGPXGPMFKCGPMKK", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "ABCGPX");
            Assert.AreEqual(last, "GPMFKCGPMKK");
        }

        [Test]
        public static void TestNTerminusProtease()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("Asp-N", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("PADDMSKDPDMMAASMDJSSM", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "PA");
            Assert.AreEqual(last, "DJSSM");
        }

        [Test]
        public static void TestWrongSyntax()
        {
            Assert.Throws<MzLibException>(() =>
            {
                var protease = DigestionMotif.ParseDigestionMotifsFromString("X[Y,P]");
                Assert.Fail("Exception shold be thrown for incorrect syntax.");
            });
        }

        [Test]
        public static void TestCutIndexDifferentSyntax()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("K|[P]"); // same as K[P]|
            Protease protease = new Protease("lys-c", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            DigestionParams myDigestionParams = new DigestionParams("lys-c", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("PROKPKMKP", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "PROKPK");
            Assert.AreEqual(last, "MKP");
        }

        [Test]
        public static void TestEndSequenceCTerm()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("chymotrypsin|P", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("AASFPWDJSSMF", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "AASFPW");
            Assert.AreEqual(last, "DJSSMF");
        }

        [Test]
        public static void TestNonSpecificProtease()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("non-specific", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("PRO", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            Assert.AreEqual(myPeptides.Count(), 3);
        }

        [Test]
        public static void TestSpecificProteaseWriting()
        {
            //check that the specific protease is the one written for indexing
            //This is needed for speedy non-specific searches to have the stratified full/semi/none peptide cleavages
            //If the protease is written instead of the specific protease, then the specific protease information is lost upon deserialization.
            //check for nonspecific
            DigestionParams dp = new DigestionParams(protease: "Arg-C", searchModeType: CleavageSpecificity.None);
            string proteaseString = dp.ToString().Split(',')[6];
            Assert.IsTrue(proteaseString.Equals("Arg-C"));
            //Check for semi
            dp = new DigestionParams(protease: "Arg-C", searchModeType: CleavageSpecificity.Semi);
            proteaseString = dp.ToString().Split(',')[6];
            Assert.IsTrue(proteaseString.Equals("Arg-C"));
            //check for normal
            dp = new DigestionParams(protease: "Arg-C"); //default searchModeType is Full
            proteaseString = dp.ToString().Split(',')[6];
            Assert.IsTrue(proteaseString.Equals("Arg-C"));
        }

        [Test]
        public static void TestEndSequenceNTerm()
        {
            var empty = new List<Modification>();
            DigestionParams myDigestionParams = new DigestionParams("Lys-N", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("KKPROTEIN", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(first, "K");
            Assert.AreEqual(myPeptides.Count(), 2);
        }

        [Test]
        public static void TestOneMotifMultiplePreventing()
        {
            var empty = new List<Modification>();
            var digestionmotifs = DigestionMotif.ParseDigestionMotifsFromString("N[M]|,N[C]|,N[A]|");
            Protease customProtease = new Protease("custom", CleavageSpecificity.Full, "", "", digestionmotifs);
            ProteaseDictionary.Dictionary.Add(customProtease.Name, customProtease);
            DigestionParams myDigestionParams = new DigestionParams("custom", minPeptideLength: 1, maxMissedCleavages: 0);
            // create a protein
            Protein myProtein = new Protein("PRONFNMMHFHAA", "myAccession");
            // digest it into peptides
            var myPeptides = myProtein.Digest(myDigestionParams, empty, empty).ToList();
            string first = myPeptides.First().ToString();
            string last = myPeptides.Last().ToString();
            Assert.AreEqual(myPeptides.Count, 2);
            Assert.AreEqual(first, "PRON");
        }

        [Test]
        public static void TestNterminalProteolysis()
        {
            //N-terminal digestion
            Protein p = new Protein("MPEPTIDE", "P12345");
            int fullProteinOneBasedBegin = 1;
            int fullProteinOneBasedEnd = 8;
            bool addNterminalDegestionTruncations = true;
            bool addCterminalDigestionTruncations = false;
            int minProductBaseSequenceLength = 2;
            int lengthOfProteolysis = 3;
            string proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            List<TruncationProduct> products = p.TruncationProducts.ToList();
            Assert.AreEqual(4, products.Count);
            List<string> productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            List<string> expectedProductSequences = new List<string> { "PEPTIDE", "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("MPEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 8;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = false;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName); products = p.TruncationProducts.ToList();
            Assert.AreEqual(4, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> {"PEPTIDE", "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            
            p = new Protein("PEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 7;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = false;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            products = p.TruncationProducts.ToList();
            Assert.AreEqual(3, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("PEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 7;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = false;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName); products = p.TruncationProducts.ToList();
            Assert.AreEqual(3, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("MPEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 8;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = false;
            minProductBaseSequenceLength = 6;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            products = p.TruncationProducts.ToList();
            Assert.AreEqual(2, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "PEPTIDE", "EPTIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
        }

        [Test]
        public static void TestCterminalProteolysis()
        {
            //N-terminal digestion
            Protein p = new Protein("MPEPTIDE", "P12345");
            int fullProteinOneBasedBegin = 1;
            int fullProteinOneBasedEnd = 8;
            bool addNterminalDegestionTruncations = false;
            bool addCterminalDigestionTruncations = true;
            int minProductBaseSequenceLength = 2;
            int lengthOfProteolysis = 3;
            string proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            List<TruncationProduct> products = p.TruncationProducts.ToList();
            Assert.AreEqual(6, products.Count);
            List<string> productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            List<string> expectedProductSequences = new List<string> { "MPEPTID", "MPEPTI", "MPEPT", "PEPTID", "PEPTI", "PEPT" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("MPEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 8;
            addNterminalDegestionTruncations = false;
            addCterminalDigestionTruncations = true;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName); products = p.TruncationProducts.ToList();
            Assert.AreEqual(6, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "MPEPTID", "MPEPTI", "MPEPT", "PEPTID", "PEPTI", "PEPT" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("PEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 7;
            addNterminalDegestionTruncations = false;
            addCterminalDigestionTruncations = true;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            products = p.TruncationProducts.ToList();
            Assert.AreEqual(3, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "PEPTID", "PEPTI", "PEPT" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("PEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 7;
            addNterminalDegestionTruncations = false;
            addCterminalDigestionTruncations = true;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            products = p.TruncationProducts.ToList();
            Assert.AreEqual(3, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "PEPTID", "PEPTI", "PEPT" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
        }

        [Test]
        public static void TestProteolysisBothTermini()
        {
            Protein p = new Protein("MPEPTIDE", "P12345");
            int fullProteinOneBasedBegin = 1;
            int fullProteinOneBasedEnd = 8;
            bool addNterminalDegestionTruncations = true;
            bool addCterminalDigestionTruncations = true;
            int minProductBaseSequenceLength = 2;
            int lengthOfProteolysis = 3;
            string proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            List<TruncationProduct> products = p.TruncationProducts.ToList();
            Assert.AreEqual(10, products.Count);
            List<string> productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            List<string> expectedProductSequences = new List<string> { "MPEPTID", "MPEPTI", "MPEPT", "PEPTID", "PEPTI", "PEPT", "PEPTIDE", "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("MPEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 8;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = true;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName); products = p.TruncationProducts.ToList();
            Assert.AreEqual(10, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "MPEPTID", "MPEPTI", "MPEPT", "PEPTIDE", "PEPTID", "PEPTI", "PEPT", "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("PEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 7;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = true;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            products = p.TruncationProducts.ToList();
            Assert.AreEqual(6, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "PEPTID", "PEPTI", "PEPT", "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("PEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 7;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = true;
            minProductBaseSequenceLength = 2;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName);
            products = p.TruncationProducts.ToList();
            Assert.AreEqual(6, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> { "PEPTID", "PEPTI", "PEPT", "EPTIDE", "PTIDE", "TIDE" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
            p = new Protein("MPEPTIDE", "P12345");
            fullProteinOneBasedBegin = 1;
            fullProteinOneBasedEnd = 8;
            addNterminalDegestionTruncations = true;
            addCterminalDigestionTruncations = true;
            minProductBaseSequenceLength = 6;
            lengthOfProteolysis = 3;
            proteolyisisProductName = "truncation";
            p.AddTruncationsToExistingProteolysisProducts(fullProteinOneBasedBegin, fullProteinOneBasedEnd, addNterminalDegestionTruncations, addCterminalDigestionTruncations, minProductBaseSequenceLength, lengthOfProteolysis, proteolyisisProductName); products = p.TruncationProducts.ToList();
            Assert.AreEqual(5, products.Count);
            productSequences = new List<string>();
            foreach (TruncationProduct product in products)
            {
                productSequences.Add(p.BaseSequence.Substring((int)product.OneBasedBeginPosition - 1, (int)product.OneBasedEndPosition - (int)product.OneBasedBeginPosition + 1));
            }
            expectedProductSequences = new List<string> {"PEPTIDE", "EPTIDE", "PEPTID", "MPEPTID", "MPEPTI" };
            CollectionAssert.AreEquivalent(expectedProductSequences, productSequences);
        }

        [Test]
        public static void TestProteoformsCleavedOnce()
        {
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "P08709.xml");
            Protein insulin = ProteinDbLoader.LoadProteinXML(xmlDatabase, true, DecoyType.None, null, false, null, out var unknownModifications)[0];
            insulin.CleaveOnceBetweenProteolysisProducts();
            List<string> productNames = insulin.TruncationProducts.Select(t => t.Type).ToList();
            Assert.AreEqual(8, productNames.Count);
            Assert.IsTrue(productNames.Contains("C-terminal Portion of Singly Cleaved Protein(21-466)"));
            Assert.IsTrue(productNames.Contains("N-terminal Portion of Singly Cleaved Protein(1-60)"));
            Assert.IsTrue(productNames.Contains("C-terminal Portion of Singly Cleaved Protein(61-466)"));
            Assert.IsTrue(productNames.Contains("N-terminal Portion of Singly Cleaved Protein(1-212)"));
        }

        [Test]
        public static void TestProteoformsCleavedOnceLong()
        {
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "P08709.xml");
            Protein insulin = ProteinDbLoader.LoadProteinXML(xmlDatabase, true, DecoyType.None, null, false, null, out var unknownModifications)[0];
            insulin.CleaveOnceBetweenProteolysisProducts(minimumProductLength: 70);
            List<string> productNames = insulin.TruncationProducts.Select(t => t.Type).ToList();
            Assert.AreEqual(7, productNames.Count);
            Assert.IsTrue(productNames.Contains("C-terminal Portion of Singly Cleaved Protein(21-466)"));
            Assert.IsTrue(!productNames.Contains("N-terminal Portion of Singly Cleaved Protein(1-60)"));
            Assert.IsTrue(productNames.Contains("C-terminal Portion of Singly Cleaved Protein(61-466)"));
            Assert.IsTrue(productNames.Contains("N-terminal Portion of Singly Cleaved Protein(1-212)"));
        }

        [Test]
        public static void TestProteolyticDigestion()
        {
            Protein humanInsulin = new Protein("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "P01308",
            proteolysisProducts: new List<TruncationProduct>
            {
                new TruncationProduct(1, 24, ""),
                new TruncationProduct(25, 54, ""),
                new TruncationProduct(57, 87, ""),
                new TruncationProduct(90, 110, "")
            });
            DigestionParams dp = new DigestionParams(maxMissedCleavages: 10, minPeptideLength: 1, maxPeptideLength: 120); //this should allow for all peptides to be generated
            List<PeptideWithSetModifications> pwsms = humanInsulin.Digest(dp, null, null).ToList();
            HashSet<PeptideWithSetModifications> hashset = new HashSet<PeptideWithSetModifications>(pwsms);
            //check that all the full length proteolysis products were made
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 1 && x.OneBasedEndResidueInProtein == 24));
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 25 && x.OneBasedEndResidueInProtein == 54));
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 57 && x.OneBasedEndResidueInProtein == 87));
            Assert.IsTrue(pwsms.Any(x => x.OneBasedStartResidueInProtein == 90 && x.OneBasedEndResidueInProtein == 110));
            //check that all the correct peptides were made
            Assert.IsTrue(hashset.Count == 52);
            //check that there are no duplicates
            Assert.IsTrue(pwsms.Count == hashset.Count);
            //Speedy semi specific test
            DigestionParams speedySemiN = new DigestionParams("trypsin", 10, 29, 30, 1024, InitiatorMethionineBehavior.Retain, 2, CleavageSpecificity.Semi, FragmentationTerminus.N);
            DigestionParams speedySemiC = new DigestionParams("trypsin", 10, 29, 30, 1024, InitiatorMethionineBehavior.Retain, 2, CleavageSpecificity.Semi, FragmentationTerminus.C);
            List<PeptideWithSetModifications> pwsmsN = humanInsulin.Digest(speedySemiN, null, null).ToList();
            List<PeptideWithSetModifications> pwsmsC = humanInsulin.Digest(speedySemiC, null, null).ToList();
            Assert.IsTrue(pwsmsN.Count == 7);
            Assert.IsTrue(pwsmsC.Count == 9);
            Assert.IsFalse(pwsmsN.Any(x => x.Length > speedySemiN.MaxLength));
            Assert.IsFalse(pwsmsC.Any(x => x.Length > speedySemiC.MaxLength));
            Assert.IsFalse(pwsmsN.Any(x => x.Length < speedySemiN.MinLength));
            Assert.IsFalse(pwsmsC.Any(x => x.Length < speedySemiC.MinLength));
        }

        [Test]
        public void TestDigestionParamsMaskedProperties()
        {
            var digestionParams = new DigestionParams();
            digestionParams.MinPeptideLength = 1;
            Assert.That(digestionParams.MinLength, Is.EqualTo(digestionParams.MinPeptideLength));

            digestionParams.MaxPeptideLength = 2;
            Assert.That(digestionParams.MaxLength, Is.EqualTo(digestionParams.MaxPeptideLength));


            digestionParams.MaxModsForPeptide = 3;
            Assert.That(digestionParams.MaxMods, Is.EqualTo(digestionParams.MaxModsForPeptide));
        }

        private class TestDigestionAgent : DigestionAgent
        {
            public TestDigestionAgent(string name, CleavageSpecificity cleavageSpecificity, List<DigestionMotif> motifList, Modification cleavageMod)
                : base(name, cleavageSpecificity, motifList, cleavageMod)
            {
            }
        }

        [Test]
        public void Equals_SameName_ReturnsTrue()
        {
            var agent1 = ProteaseDictionary.Dictionary["trypsin"];
            var agent2 = ProteaseDictionary.Dictionary["trypsin"];

            Assert.That(agent1.Equals(agent2), Is.True);
        }

        [Test]
        public void Equals_DifferentName_ReturnsFalse()
        {
            var agent1 = ProteaseDictionary.Dictionary["trypsin"];
            var agent2 = ProteaseDictionary.Dictionary["Arg-C"];

            Assert.That(agent1.Equals(agent2), Is.False);
        }

        [Test]
        public void GetHashCode_SameName_ReturnsSameHashCode()
        {
            var agent1 = ProteaseDictionary.Dictionary["trypsin"];
            var agent2 = ProteaseDictionary.Dictionary["trypsin"];

            Assert.That(agent1.GetHashCode(), Is.EqualTo(agent2.GetHashCode()));
        }

        [Test]
        public void GetHashCode_DifferentName_ReturnsDifferentHashCode()
        {
            var agent1 = ProteaseDictionary.Dictionary["trypsin"];
            var agent2 = ProteaseDictionary.Dictionary["Arg-C"];

            Assert.That(agent1.GetHashCode(), Is.Not.EqualTo(agent2.GetHashCode()));
        }
    }
}