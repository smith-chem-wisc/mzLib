using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.InternalFileReading
{
    /// <summary>
    /// MetaMorpheus's DIA precursor table (AllDiaPrecursors.tsv): one row per precursor per run. The writer and this reader
    /// live together so the columns cannot drift (dataRepo's DATAREPO-65). The schema is the one agreed with dataRepo
    /// (dia 003): the run key verbatim, labels T/D/ET/ED/C, level-named q-values, minutes and iRT in their own columns,
    /// an apex scan number, and an empty quantity when none is reported.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestDiaPrecursorFile
    {
        private static string Folder => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestDiaPrecursorFile");

        [SetUp]
        public void SetUp() => Directory.CreateDirectory(Folder);

        [TearDown]
        public void TearDown() => Directory.Delete(Folder, true);

        private static DiaPrecursorFromTsv Row(string run, string label, double? quantity) => new()
        {
            FileName = run,
            FullSequence = "PEPM[Common Variable:Oxidation on M]TIDE",
            BaseSequence = "PEPMTIDE",
            PrecursorCharge = 2,
            PrecursorMz = 466.7012,
            ProteinAccession = "P12345|Q67890",
            Label = label,
            Score = 0.9871,
            QValuePrecursorRun = 0.0012,
            QValuePrecursorGlobal = 0.0021,
            LibraryIrt = 41.5,
            ApexIrt = 42.25,
            ApexRtMin = 17.125,
            ApexScanNumber = 30411,
            PrecursorQuantity = quantity,
        };

        /// <summary>
        /// Every column survives a write and a read, including a PRIDE-style run name whose extensions would be stripped
        /// by the psmtsv name handling, an unreported quantity (empty, not 0) and whether contaminants were assessed.
        /// </summary>
        [Test]
        public void EveryColumnRoundTrips()
        {
            var rows = new List<DiaPrecursorFromTsv>
            {
                Row("X.raw.thermo.raw", "T", 12345.5),
                Row("X.raw.thermo.raw", "D", null),
                Row("run 2.mzML", "ET", 0),
                Row("run 2.mzML", "C", 7.25),
            };
            string path = Path.Combine(Folder, "AllDiaPrecursors.tsv");
            new DiaPrecursorFile(path, rows, contaminantsAssessed: false).WriteResults(path);

            var read = new DiaPrecursorFile(path);
            read.LoadResults();

            Assert.That(read.ContaminantsAssessed, Is.False);
            Assert.That(read.Results.Count, Is.EqualTo(rows.Count));
            for (int i = 0; i < rows.Count; i++)
                Assert.That(read.Results[i], Is.EqualTo(rows[i]), $"row {i}");
            Assert.That(read.Results[1].PrecursorQuantity, Is.Null, "an unreported quantity stays unreported");
            Assert.That(read.Results[2].PrecursorQuantity, Is.EqualTo(0), "a reported 0 stays 0");
        }

        /// <summary>
        /// The file says whether a contaminant database was searched on its first line, so "no contaminants" is never
        /// true by default; the header names each q-value's level and the unit of each retention column.
        /// </summary>
        [Test]
        public void TheFileStatesContaminantsAndUnits()
        {
            string path = Path.Combine(Folder, "AllDiaPrecursors.tsv");
            new DiaPrecursorFile(path, [Row("a.raw", "T", 1)], contaminantsAssessed: true).WriteResults(path);
            var lines = File.ReadAllLines(path);

            Assert.That(lines[0], Is.EqualTo("# contaminants: assessed"));
            var header = lines[1].Split('\t');
            Assert.That(header, Is.EqualTo(new[]
            {
                "File Name", "Full Sequence", "Base Sequence", "Precursor Charge", "Precursor MZ", "Protein Accession",
                "Decoy/Contaminant/Target", "Score", "QValue_Precursor_Run", "QValue_Precursor_Global", "LibraryIrt", "ApexIrt",
                "ApexRtMin", "ApexScanNumber", "Precursor Quantity",
            }));
            Assert.That(lines[2].Split('\t')[0], Is.EqualTo("a.raw"));

            var read = new DiaPrecursorFile(path);
            read.LoadResults();
            Assert.That(read.ContaminantsAssessed, Is.True);
        }

        /// <summary>The file is recognised by its name, as MetaMorpheus writes it, and opens through the reader factory.</summary>
        [Test]
        public void TheFileTypeIsRecognised()
        {
            string path = Path.Combine(Folder, "AllDiaPrecursors.tsv");
            new DiaPrecursorFile(path, [Row("a.raw", "T", 1)], contaminantsAssessed: false).WriteResults(path);

            Assert.That(path.ParseFileType(), Is.EqualTo(SupportedFileType.DiaPrecursorTsv));
            var file = FileReader.ReadFile<DiaPrecursorFile>(path);
            Assert.That(file.Results.Single().FileName, Is.EqualTo("a.raw"));
        }

        /// <summary>A file without the contaminant line is refused rather than read as "not assessed".</summary>
        [Test]
        public void AFileWithoutTheContaminantLineIsRefused()
        {
            string path = Path.Combine(Folder, "AllDiaPrecursors.tsv");
            new DiaPrecursorFile(path, [Row("a.raw", "T", 1)], contaminantsAssessed: false).WriteResults(path);
            File.WriteAllLines(path, File.ReadAllLines(path).Skip(1));

            Assert.Throws<MzLibUtil.MzLibException>(() => new DiaPrecursorFile(path).LoadResults());
        }
    }
}
