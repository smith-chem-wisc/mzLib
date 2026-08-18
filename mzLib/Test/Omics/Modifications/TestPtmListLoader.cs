using MzLibUtil;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Modifications;
using ModificationLoader = global::Omics.Modifications.IO.ModificationLoader;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.Omics.Modifications
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestPtmListLoader
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
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
        public static void SampleModFileLoading()
        {
            PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "sampleModFile.txt"), out var errors);
        }

        [Test]
        [TestCase("CommonArtifacts.txt", 33)]
        [TestCase("CommonBiological.txt", 35)]
        public static void Test_ReadAllModsFromFile(string filename, int modCount)
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", filename);
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation, out var errors);
            Assert.AreEqual(modCount, a.Count());
        }

        [Test]
        [TestCase("CommonArtifacts.txt")]
        [TestCase("CommonBiological.txt")]
        public static void Test_ModsFromFileAreSorted(string filename)
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", filename);
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation, out var errors);

            string id1 = a.First().IdWithMotif;
            foreach (string modId in a.Select(m => m.IdWithMotif))
            {
                Assert.GreaterOrEqual(modId, id1);
            }
        }

        [Test]
        public static void Test_ModsWithComments()
        {
            string testModificationsFileLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "ModificationTests", @"ModsWithComments.txt");
            var a = PtmListLoader.ReadModsFromFile(testModificationsFileLocation, out var errors).ToList();
            Assert.AreEqual(4, a.Select(m => m.IdWithMotif).ToList().Count);

            Assert.AreEqual("Deamidation on N", a[0].IdWithMotif.ToString());
            Assert.AreEqual("Sodium on D", a[2].IdWithMotif.ToString());//this has trailing whitespace that shouldn't be in the name

            //Make sure comments are okay on DR key and that key value pairs are still split correctly
            Modification someMod = a[2];
            Modification test = someMod;
            var residValueTest = test.DatabaseReference.First().Value.First();
            var residKeyTest = test.DatabaseReference.First().Key;
            Assert.AreEqual("RESID", residKeyTest);
            Assert.AreEqual("AA0441", residValueTest);
        }

        [Test]
        [TestCase("sampleModFileFail1.txt")] // TG is not valid
        [TestCase("sampleModFileFail2.txt")]
        [TestCase("sampleModFileFail5.txt")] // ID is missing
        [TestCase("sampleModFileFail6.txt")] // modification type is missing
        [TestCase("sampleModFileFail_missingPosition.txt")] // missing position
        [TestCase("sampleModFileFail_missingChemicalFormulaAndMonoisotopicMass.txt")]
        public static void SampleModFileLoadingFail1General(string filename)
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), out var errors).ToList();
            Assert.AreEqual(0, a.Count);
        }

        [Test]
        public static void PTMListLoader_ModWithComments_Equals_ModWithoutComments()
        {
            var a = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_Comments.txt"), out var errors).ToList();
            var b = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SampleMod_NoComments.txt"), out var errors2).ToList();
            Assert.IsTrue(a.First().Equals(b.First()));
        }

        [Test]
        [TestCase("sampleModFileFail3.txt", "Input string for chemical formula was in an incorrect format: $%&$%")]
        [TestCase("m.txt", "0 or 238.229666 is not a valid monoisotopic mass")]
        public static void SampleModFileLoadingFail3General(string filename, string errorMessage)
        {
            NUnit.Framework.Assert.That(() =>
                PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), out var errors).ToList(),
                Throws.TypeOf<MzLibException>().With.Property("Message").EqualTo(errorMessage));
        }

        [Test]
        [TestCase("sampleModFileDouble.txt")]
        [TestCase("sampleModFileDouble2.txt")]
        public static void CompactFormReadingGeneral(string filename)
        {
            Assert.AreEqual(2, PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), out var errors).Count());
        }

        [Test]
        public static void TestReadingInvalidModifications()
        {
            string modText = "ID   mod\r\nMT   type\r\nPP   Anywhere.\r\nTG   X\r\nCF   H1\r\n" + @"//";
            var mods = PtmListLoader.ReadModsFromString(modText, out var errors1);
            Assert.AreEqual(1, mods.Count());
            Assert.AreEqual(0, errors1.Count);

            modText = "ID   mod\r\nMT   type\r\nPP   Anywhere.\r\nTG   X\r\n" + @"//"; // no monoisotopic mass, so invalid
            mods = PtmListLoader.ReadModsFromString(modText, out var errors2);
            Assert.AreEqual(0, mods.Count());
            Assert.AreEqual(1, errors2.Count);
            Assert.IsFalse(errors2.Single().Item1.ValidModification);
            Assert.IsTrue(errors2.Single().Item2.Split(new[] { "\r\n" }, StringSplitOptions.None).Any(x => x.StartsWith("#"))); // has an error comment
        }

        [Test]
        public static void TestReadingIdWithMotif()
        {
            string modText = "ID   Detached EVK or XleDK\r\nPP   Peptide N-terminal.\r\nTG   evkX or vekX or ldkX or dlkX or idkX or dikX\r\nMT   Detached\r\nNL   C16H28N4O5\r\nCF   C16H28N4O5\r\n" + @"//";

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "detacher.txt");
            File.WriteAllLines(path, new string[] { modText });

            var mods = PtmListLoader.ReadModsFromFile(path, out var errors).ToList();
            var motifs = mods.Select(p => p.Target.ToString()).Distinct().ToList();
            var ids = mods.Select(p => p.IdWithMotif).Distinct().ToList();

            NUnit.Framework.Assert.That(mods.Count == 6);
            NUnit.Framework.Assert.That(motifs.Count == 6);
            NUnit.Framework.Assert.That(ids.Count == 6);
        }

        [Test]
        public static void TestInvalidModTypeError()
        {
            string mod = "ID   Deamidation\r\nTG   N or Q\r\nPP   Anywhere.\r\nMT   Mod:\r\nCF   H-1 N-1 O1\r\n//";
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestInvalidModTypeError\ptmlist.txt");
            Directory.CreateDirectory(Directory.GetParent(filepath).FullName);
            File.WriteAllLines(filepath, new string[] { mod });

            var ptms = PtmListLoader.ReadModsFromFile(filepath, out var warnings).ToList();
            NUnit.Framework.Assert.That(ptms.Count == 0);
            NUnit.Framework.Assert.That(warnings.Count == 2);

            NUnit.Framework.Assert.That(warnings.First().Item2.Contains("Modification type cannot contain ':'!"));
        }

        // Issue #353: a field code with no value used to be dereferenced without a null check. Existing tests
        // only covered an absent line, which is a different path. Both loaders are covered because
        // PtmListLoader duplicates the parsing loop rather than delegating to ModificationLoader.
        private const string ValidMod = "ID   testmod\r\nMT   type\r\nPP   Anywhere.\r\nTG   X\r\nCF   H1\r\n//";

        [Test]
        // padded to the usual column and bare, because the loader trims before it slices the value off
        [TestCase("ID")] [TestCase("MT")] [TestCase("PP")] [TestCase("TG")] [TestCase("CF")]
        [TestCase("MM")] [TestCase("DR")] [TestCase("TR")] [TestCase("KW")] [TestCase("NL")]
        [TestCase("DI")] [TestCase("FT")] [TestCase("AC")]
        public static void TestFieldCodeWithNoValueDoesNotThrow(string fieldCode)
        {
            foreach (string emptyLine in new[] { fieldCode, fieldCode + "   " })
            {
                string modText = ValidMod.Replace("//", emptyLine + "\r\n//");

                foreach ((string name, ReadModsFromString read) in Loaders)
                {
                    NUnit.Framework.Assert.DoesNotThrow(() => read(modText, out _).ToList(),
                        $"{name} threw on an empty '{fieldCode}' value");
                }
            }
        }

        // Both loaders, every time: PtmListLoader duplicates the parsing loop rather than delegating to
        // ModificationLoader, so a guard added to one is not a guard on the other.
        private delegate IEnumerable<Modification> ReadModsFromString(string text, out List<(Modification, string)> errors);

        private static IEnumerable<(string Name, ReadModsFromString Read)> Loaders =>
            new (string, ReadModsFromString)[]
            {
                ("ModificationLoader", ModificationLoader.ReadModsFromString),
#pragma warning disable CS0618 // obsolete, still public, still crashes identically
                ("PtmListLoader", PtmListLoader.ReadModsFromString),
#pragma warning restore CS0618
            };

        private static Dictionary<string, IList<string>> ReferenceFor(Modification mod, string fieldCode) =>
            fieldCode == "DR" ? mod.DatabaseReference : mod.TaxonomicRange;

        // DR and TR index straight into the split result, so a value that is not a "key; value" pair threw.
        // Skipping the line is only correct if it skips nothing else: the record and its other fields still
        // have to load, and a well-formed pair still has to parse. Without that second half the guard could
        // be discarding the field wholesale and the first half would not notice.
        [Test]
        [TestCase("DR")]
        [TestCase("TR")]
        public static void TestReferenceWithoutAPairIsSkippedAndAValidPairStillParses(string fieldCode)
        {
            string malformed = ValidMod.Replace("//", fieldCode + "   onlykey\r\n//");
            string wellFormed = ValidMod.Replace("//", fieldCode + "   somekey; someval\r\n//");

            foreach ((string name, ReadModsFromString read) in Loaders)
            {
                var fromMalformed = read(malformed, out var malformedErrors).ToList();

                NUnit.Framework.Assert.That(fromMalformed.Count, Is.EqualTo(1),
                    $"{name} discarded the whole record over a malformed {fieldCode}");
                NUnit.Framework.Assert.That(malformedErrors.Count, Is.EqualTo(0));
                NUnit.Framework.Assert.That(fromMalformed[0].IdWithMotif, Is.EqualTo("testmod on X"));
                NUnit.Framework.Assert.That(fromMalformed[0].MonoisotopicMass, Is.Not.Null,
                    $"{name} lost a field that came after the malformed {fieldCode}");
                NUnit.Framework.Assert.That(ReferenceFor(fromMalformed[0], fieldCode), Is.Null,
                    $"{name} half-parsed a malformed {fieldCode} into the record");

                var fromWellFormed = read(wellFormed, out var wellFormedErrors).ToList();

                NUnit.Framework.Assert.That(fromWellFormed.Count, Is.EqualTo(1));
                NUnit.Framework.Assert.That(wellFormedErrors.Count, Is.EqualTo(0));
                NUnit.Framework.Assert.That(ReferenceFor(fromWellFormed[0], fieldCode)["somekey"],
                    Is.EquivalentTo(new[] { "someval" }),
                    $"{name} stopped parsing a well-formed {fieldCode} pair");
            }
        }

        // The guard's contract is that an empty value is treated as an absent line. So a blank mandatory
        // field has to fail the same validation a missing one fails, with the same message, rather than
        // becoming a new error mode. Compared against the absent-line result rather than a hardcoded
        // string, so this stays true if the message is ever reworded.
        [Test]
        [TestCase("TG", "ID   m|MT   t|PP   Anywhere.|TG   X|CF   H1")]
        [TestCase("CF", "ID   m|MT   t|PP   Anywhere.|TG   X|CF   H1")]
        [TestCase("MM", "ID   m|MT   t|PP   Anywhere.|TG   X|MM   57.02146")]
        public static void TestBlankMandatoryFieldIsTreatedAsAnAbsentLine(string fieldCode, string record)
        {
            string[] lines = record.Split('|');
            string blank = string.Join("\r\n", lines.Select(x => x.StartsWith(fieldCode) ? fieldCode + "   " : x)) + "\r\n//";
            string absent = string.Join("\r\n", lines.Where(x => !x.StartsWith(fieldCode))) + "\r\n//";

            foreach ((string name, ReadModsFromString read) in Loaders)
            {
                var fromBlank = read(blank, out var blankErrors).ToList();
                var fromAbsent = read(absent, out var absentErrors).ToList();

                NUnit.Framework.Assert.That(fromBlank.Count, Is.EqualTo(0),
                    $"{name} accepted a record whose mandatory {fieldCode} was blank");
                NUnit.Framework.Assert.That(fromAbsent.Count, Is.EqualTo(0));
                NUnit.Framework.Assert.That(blankErrors.Select(x => x.Item2), Is.EqualTo(absentErrors.Select(x => x.Item2)),
                    $"{name} reported a blank {fieldCode} differently from an absent {fieldCode}");
            }
        }

        // The reachable form of the bug: these are values mzLib's own writer emits, so a database it wrote
        // could not be read back. An empty keyword writes "KW   " and an empty reference writes "DR   ; ".
        // The empty entry itself is dropped rather than preserved - the assertion is that the file loads at
        // all, not that the round trip is lossless.
        [Test]
        public static void TestModificationsMzLibWritesCanBeReadBack()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);

            var emptyKeyword = new Modification("m", null, "mt", null, motif, "Anywhere.", null, 10,
                _keywords: new List<string> { "" });
            var emptyReference = new Modification("m", null, "mt", null, motif, "Anywhere.", null, 10,
                _databaseReference: new Dictionary<string, IList<string>> { { "", new List<string> { "" } } });

            foreach (var modification in new[] { emptyKeyword, emptyReference })
            {
                string written = modification.ToString() + "\r\n//";

                foreach ((string name, ReadModsFromString read) in Loaders)
                {
                    var reread = read(written, out var errors).ToList();

                    NUnit.Framework.Assert.That(reread.Count, Is.EqualTo(1),
                        $"{name} could not read back a modification mzLib wrote:\r\n{written}");
                    NUnit.Framework.Assert.That(errors.Count, Is.EqualTo(0));
                    NUnit.Framework.Assert.That(reread[0].ValidModification, Is.True);
                    NUnit.Framework.Assert.That(reread[0].IdWithMotif, Is.EqualTo("m on X"));
                }
            }
        }
    }
}