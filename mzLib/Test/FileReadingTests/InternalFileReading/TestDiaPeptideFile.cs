using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.InternalFileReading
{
    /// <summary>
    /// MetaMorpheus's DIA peptide table (AllDiaPeptides.tsv): one row per peptide across the runs searched together, with
    /// the peptide-level q-value. It follows the precursor table's conventions (dataRepo's eight points, dia 004): labels
    /// T/D/ET/ED/C, a level-named q-value, the run key verbatim, and a first line saying whether contaminants were assessed.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestDiaPeptideFile
    {
        private static string Folder => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestDiaPeptideFile");

        [SetUp]
        public void SetUp() => Directory.CreateDirectory(Folder);

        [TearDown]
        public void TearDown() => Directory.Delete(Folder, true);

        private static DiaPeptideFromTsv Row(string run, string label) => new()
        {
            FullSequence = "PEPM[Common Variable:Oxidation on M]TIDE",
            BaseSequence = "PEPMTIDE",
            ProteinAccession = "P12345|Q67890",
            Label = label,
            Score = 0.9871,
            QValuePeptideGlobal = 0.0021,
            FileName = run,
            PrecursorCharge = 3,
        };

        /// <summary>
        /// Every column survives a write and a read, including a PRIDE-style run name whose extensions the psmtsv name
        /// handling would strip, and whether contaminants were assessed.
        /// </summary>
        [Test]
        public void EveryColumnRoundTrips()
        {
            var rows = new List<DiaPeptideFromTsv>
            {
                Row("X.raw.thermo.raw", "T"),
                Row("X.raw.thermo.raw", "D"),
                Row("run 2.mzML", "ET"),
                Row("run 2.mzML", "C"),
            };
            string path = Path.Combine(Folder, "AllDiaPeptides.tsv");
            new DiaPeptideFile(path, rows, contaminantsAssessed: true).WriteResults(path);

            var read = new DiaPeptideFile(path);
            read.LoadResults();

            Assert.That(read.ContaminantsAssessed, Is.True);
            Assert.That(read.Results, Is.EqualTo(rows));
        }

        /// <summary>
        /// The first line states whether a contaminant database was searched; the header names the q-value's level, and
        /// the best precursor's run and charge come last, as attributes of the peptide rather than its key.
        /// </summary>
        [Test]
        public void TheFileStatesContaminantsAndTheQValueLevel()
        {
            string path = Path.Combine(Folder, "AllDiaPeptides.tsv");
            new DiaPeptideFile(path, [Row("a.raw", "T")], contaminantsAssessed: false).WriteResults(path);
            var lines = File.ReadAllLines(path);

            Assert.That(lines[0], Is.EqualTo("# contaminants: not assessed"));
            Assert.That(lines[1].Split('\t'), Is.EqualTo(new[]
            {
                "Full Sequence", "Base Sequence", "Protein Accession", "Decoy/Contaminant/Target", "Score",
                "QValue_Peptide_Global", "File Name", "Precursor Charge",
            }));

            var read = new DiaPeptideFile(path);
            read.LoadResults();
            Assert.That(read.ContaminantsAssessed, Is.False);
        }

        /// <summary>A file without the contaminant line is refused rather than read as "not assessed".</summary>
        [Test]
        public void AFileWithoutTheContaminantLineIsRefused()
        {
            string path = Path.Combine(Folder, "AllDiaPeptides.tsv");
            new DiaPeptideFile(path, [Row("a.raw", "T")], contaminantsAssessed: false).WriteResults(path);
            File.WriteAllLines(path, File.ReadAllLines(path).Skip(1));

            Assert.Throws<MzLibUtil.MzLibException>(() => new DiaPeptideFile(path).LoadResults());
        }

        /// <summary>The file is recognised by the name MetaMorpheus gives it, and opens through the reader factory.</summary>
        [Test]
        public void TheFileTypeIsRecognised()
        {
            string path = Path.Combine(Folder, "AllDiaPeptides.tsv");
            new DiaPeptideFile(path, [Row("a.raw", "T")], contaminantsAssessed: false).WriteResults(path);

            Assert.That(path.ParseFileType(), Is.EqualTo(SupportedFileType.DiaPeptideTsv));
            Assert.That(FileReader.ReadFile<DiaPeptideFile>(path).Results.Single().FileName, Is.EqualTo("a.raw"));
        }

        /// <summary>A path without the table's name gets it appended, so the written file is still recognised.</summary>
        [Test]
        public void APathWithoutTheTableNameGetsItAppended()
        {
            string prefix = Path.Combine(Folder, "run1_");
            new DiaPeptideFile(prefix, [Row("a.raw", "T")], contaminantsAssessed: false).WriteResults(prefix);

            string written = prefix + "AllDiaPeptides.tsv";
            Assert.That(File.Exists(written));
            Assert.That(FileReader.ReadFile<DiaPeptideFile>(written).Results.Single().FileName, Is.EqualTo("a.raw"));
        }

        /// <summary>
        /// Columns are read by name, not by position: a table whose columns arrive in another order reads the same rows.
        /// </summary>
        [Test]
        public void ColumnsAreReadByNameInAnyOrder()
        {
            var rows = new List<DiaPeptideFromTsv> { Row("a.raw", "T"), Row("b.raw", "D") };
            string path = Path.Combine(Folder, "AllDiaPeptides.tsv");
            new DiaPeptideFile(path, rows, contaminantsAssessed: true).WriteResults(path);
            var lines = File.ReadAllLines(path);
            File.WriteAllLines(path, lines.Take(1).Concat(lines.Skip(1).Select(l => string.Join('\t', l.Split('\t').Reverse()))));

            var read = new DiaPeptideFile(path);
            read.LoadResults();

            Assert.That(File.ReadAllLines(path)[1], Does.StartWith("Precursor Charge\t"), "the columns really were reordered");
            Assert.That(read.Results, Is.EqualTo(rows));
        }

        /// <summary>
        /// The DIA protein-group table is MetaMorpheus's protein-group format, so it needs no reader of its own: its name
        /// resolves to the existing protein-group reader.
        /// </summary>
        [Test]
        public void TheDiaProteinGroupTableOpensWithTheProteinGroupReader()
        {
            Assert.That("AllDiaProteinGroups.tsv".ParseFileType(), Is.EqualTo(SupportedFileType.MetaMorpheusQuantifiedProteinGroups));
            Assert.That(SupportedFileType.MetaMorpheusQuantifiedProteinGroups.GetResultFileType(), Is.EqualTo(typeof(ProteinGroupFromTsvFile)));
        }
    }
}
