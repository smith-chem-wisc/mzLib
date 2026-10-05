using System;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.Omics.SampleInfo
{
    /// <summary>
    /// Tests for SampleExperimentalDesign, the general-purpose IExperimentalDesign implementation.
    /// The behaviour that matters to the quantification engine is the key form -- file name with
    /// extension, matched case-insensitively -- and the per-file sample order, which is what aligns
    /// samples with ISpectralMatch.Intensities.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SampleExperimentalDesignTests
    {
        private static SpectraFileInfo File(string path, string condition = "Control", int biorep = 1)
            => new SpectraFileInfo(path, condition, biorep, 1, 0);

        private static IsobaricQuantSampleInfo Channel(string path, string label, double mz, bool isReference = false)
            => new IsobaricQuantSampleInfo(path, "Control", 1, 1, 0, 1, label, mz, isReference);

        [Test]
        public void Add_KeysByFileNameWithExtension_NotByFullPath()
        {
            var design = new SampleExperimentalDesign();
            design.Add(@"C:\Data\Experiment\run1.raw", File(@"C:\Data\Experiment\run1.raw"));

            Assert.Multiple(() =>
            {
                // The engine looks up Path.GetFileName(spectralMatch.FullFilePath), so the key has to be
                // the bare name even when the design was built from full paths.
                Assert.That(design.FileNameSampleInfoDictionary.ContainsKey("run1.raw"), Is.True);
                Assert.That(design.FileNameSampleInfoDictionary.ContainsKey(@"C:\Data\Experiment\run1.raw"), Is.False);
            });
        }

        [Test]
        public void Add_AcceptsABareFileName()
        {
            var design = new SampleExperimentalDesign();
            design.Add("run1.raw", File(@"C:\Data\run1.raw"));

            Assert.That(design.FileNameSampleInfoDictionary.ContainsKey("run1.raw"), Is.True);
        }

        [Test]
        public void Lookup_IsCaseInsensitive()
        {
            var design = new SampleExperimentalDesign();
            design.Add("Sample1.raw", File(@"C:\Data\Sample1.raw"));

            // A design that names Sample1.raw should resolve data at sample1.raw. Under the default
            // ordinal comparer this lookup misses and the engine throws KeyNotFoundException.
            Assert.That(design.FileNameSampleInfoDictionary.TryGetValue("sample1.raw", out var samples), Is.True);
            Assert.That(samples, Has.Length.EqualTo(1));
        }

        [Test]
        public void Add_RejectsACaseOnlyDuplicate()
        {
            var design = new SampleExperimentalDesign();
            design.Add("Sample1.raw", File(@"C:\Data\Sample1.raw"));

            // Rejected rather than silently keeping one of the two.
            var ex = Assert.Throws<ArgumentException>(() => design.Add("sample1.raw", File(@"C:\Data\sample1.raw")));
            Assert.That(ex.Message, Does.Contain("already in this design"));
        }

        [Test]
        public void Add_PreservesChannelOrder()
        {
            const string path = @"C:\Data\tmt.raw";
            var design = new SampleExperimentalDesign();
            design.Add(path,
                Channel(path, "126", 126.12776),
                Channel(path, "127N", 127.12476),
                Channel(path, "127C", 127.13108));

            var channels = design.FileNameSampleInfoDictionary["tmt.raw"];

            // Intensities map to samples by position, so this order is the contract.
            Assert.That(channels.Cast<IsobaricQuantSampleInfo>().Select(c => c.ChannelLabel),
                Is.EqualTo(new[] { "126", "127N", "127C" }));
        }

        /// <summary>
        /// One channel of one file listed twice is refused, naming the channel and both samples.
        ///
        /// The shape is PXD040455's, a public TMT × SILAC SDRF that lists a light and a heavy sample
        /// under one reporter channel 551 times. Quantification indexes columns by sample, and a sample
        /// name takes no part in equality, so the two would merge into one column carrying one of the
        /// names over both. Refusing here is what makes leaving the name out of equality safe.
        /// </summary>
        [Test]
        public void Add_RejectsOneChannelListedTwice_NamingBothSamples()
        {
            const string path = @"C:\Data\Chip1_F2.raw";
            var light = new IsobaricQuantSampleInfo(path, "Control", 1, 1, 0, 1, "127N", 127.12476, false) { SampleName = "Chip1_F2_TMT127N_light" };
            var heavy = new IsobaricQuantSampleInfo(path, "Control", 1, 1, 0, 1, "127N", 127.12476, false) { SampleName = "Chip1_F2_TMT127N_heavy" };

            var ex = Assert.Throws<ArgumentException>(() => new SampleExperimentalDesign().Add(path, light, heavy));

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("channel 127N"));
                Assert.That(ex.Message, Does.Contain("'Chip1_F2_TMT127N_light'"));
                Assert.That(ex.Message, Does.Contain("'Chip1_F2_TMT127N_heavy'"));
            });
        }

        /// <summary>
        /// FromSamples — the entry point a design projected from SDRF will use — refuses it too.
        /// </summary>
        [Test]
        public void FromSamples_RejectsOneChannelListedTwice()
        {
            const string path = @"C:\Data\tmt.raw";
            var samples = new ISampleInfo[] { Channel(path, "126", 126.12776), Channel(path, "126", 126.12776) };

            var ex = Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.FromSamples(samples));
            Assert.That(ex.Message, Does.Contain("channel 126 of 'tmt.raw' is listed 2 times"));
        }

        /// <summary>
        /// The rule is not isobaric-only: one label-free sample listed twice for its file is refused
        /// the same way, by Add and by FromSamples, and named by its file. Master accepted both, and the
        /// second entry took the first one's column.
        /// </summary>
        [Test]
        public void Add_AndFromSamples_RejectOneLabelFreeSampleListedTwice()
        {
            const string path = @"C:\Data\run7.raw";

            var fromAdd = Assert.Throws<ArgumentException>(() => new SampleExperimentalDesign().Add(path, File(path), File(path)));
            var fromSamples = Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.FromSamples(new ISampleInfo[] { File(path), File(path) }));

            Assert.Multiple(() =>
            {
                Assert.That(fromAdd.Message, Does.Contain("sample 'run7.raw' is listed 2 times"));
                Assert.That(fromAdd.ParamName, Is.EqualTo("samples"));
                Assert.That(fromSamples.Message, Does.Contain("sample 'run7.raw' is listed 2 times"));
            });
        }

        /// <summary>
        /// An input that breaks an older rule as well as the repeat rule keeps the older rule's message.
        /// The repeat check runs last in Add, so a file added twice is still reported as already in the
        /// design, against the file parameter, even when its second samples also repeat.
        /// </summary>
        [Test]
        public void Add_ReportsAnAlreadyAddedFileBeforeARepeatedSample()
        {
            const string path = @"C:\Data\tmt.raw";
            var design = new SampleExperimentalDesign();
            design.Add(path, Channel(path, "126", 126.12776));

            var ex = Assert.Throws<ArgumentException>(
                () => design.Add(path, Channel(path, "126", 126.12776), Channel(path, "126", 126.12776)));

            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("already in this design"));
                Assert.That(ex.ParamName, Is.EqualTo("fileNameOrPath"));
            });
        }

        /// <summary>
        /// Null entries are skipped rather than counted as one sample repeated. Add refuses nulls before
        /// asking, but the engine asks the rule of any design, whose arrays may hold them.
        /// </summary>
        [Test]
        public void DescribeRepeatedSample_IgnoresNullEntries()
        {
            Assert.That(SampleExperimentalDesign.DescribeRepeatedSample(new ISampleInfo[] { null, Channel(@"C:\Data\tmt.raw", "126", 126.12776), null }),
                Is.Null);
        }

        /// <summary>
        /// The repeated channel's names are listed in ordinal order, not in the order they arrive, so the
        /// message is the same whatever order the caller's collection enumerates in. The engine asks the
        /// rule of a Dictionary's values, whose order is not guaranteed.
        /// </summary>
        [Test]
        public void DescribeRepeatedSample_ListsTheNamesInOrdinalOrder_NotArrivalOrder()
        {
            var samples = new ISampleInfo[]
            {
                new IsobaricQuantSampleInfo(@"C:\Data\tmt.raw", "Control", 1, 1, 0, 1, "126", 126.12776, false) { SampleName = "Pt3" },
                new IsobaricQuantSampleInfo(@"C:\Data\tmt.raw", "Control", 2, 1, 0, 1, "126", 126.12776, false) { SampleName = "Pt1" }
            };

            Assert.That(SampleExperimentalDesign.DescribeRepeatedSample(samples),
                Is.EqualTo("channel 126 of 'tmt.raw' is listed 2 times, as 'Pt1' and 'Pt3'"));
        }

        /// <summary>
        /// One sample on the same channel of every plex is a bridge design, not a repeat. PXD008841 puts
        /// a sample named "pool" in 131N of every plex; the channels are in different files, so they are
        /// different samples however they are named.
        ///
        /// Asked of the rule directly, not through FromSamples. FromSamples groups by file before the
        /// rule ever runs, so two plexes' pool channels never meet there and a rule that ignored the
        /// file entirely would still pass that way. The engine does hand the rule every file's samples
        /// at once, and its fixture repeats each channel label in both files, so every successful engine
        /// run is the same case end to end.
        /// </summary>
        [Test]
        public void DescribeRepeatedSample_OneSampleOnTheSameChannelOfEveryPlex_IsNotARepeat()
        {
            var samples = new ISampleInfo[]
            {
                new IsobaricQuantSampleInfo(@"C:\Data\TMTpool1_fr01.raw", "Pool", 1, 1, 1, 1, "131N", 131.13, true) { SampleName = "pool" },
                new IsobaricQuantSampleInfo(@"C:\Data\TMTpool2_fr01.raw", "Pool", 1, 1, 1, 2, "131N", 131.13, true) { SampleName = "pool" }
            };

            Assert.That(SampleExperimentalDesign.DescribeRepeatedSample(samples), Is.Null);
        }

        [Test]
        public void Add_RejectsAnEmptySampleArray()
        {
            var design = new SampleExperimentalDesign();

            var ex = Assert.Throws<ArgumentException>(() => design.Add("run1.raw"));
            Assert.That(ex.Message, Does.Contain("no samples"));
        }

        [Test]
        public void Add_RejectsANullSample()
        {
            var design = new SampleExperimentalDesign();

            // A missing channel has to be described, not omitted -- omitting shifts every later channel.
            var ex = Assert.Throws<ArgumentException>(
                () => design.Add("run1.raw", File(@"C:\Data\run1.raw"), null));
            Assert.That(ex.Message, Does.Contain("null sample"));
        }

        [Test]
        public void Add_RejectsAnEmptyFileName()
        {
            var design = new SampleExperimentalDesign();

            Assert.Throws<ArgumentException>(() => design.Add("  ", File(@"C:\Data\run1.raw")));
        }

        [Test]
        public void LabelFree_GivesOneSamplePerFile()
        {
            var files = new[]
            {
                File(@"C:\Data\a.raw", "Control"),
                File(@"C:\Data\b.raw", "Treated"),
            };

            var design = SampleExperimentalDesign.LabelFree(files);

            Assert.Multiple(() =>
            {
                Assert.That(design.FileNameSampleInfoDictionary, Has.Count.EqualTo(2));
                Assert.That(design.FileNameSampleInfoDictionary["a.raw"], Has.Length.EqualTo(1));
                Assert.That(design.FileNameSampleInfoDictionary["b.raw"].Single().Condition, Is.EqualTo("Treated"));
            });
        }

        [Test]
        public void LabelFree_RejectsARepeatedFile()
        {
            var files = new[]
            {
                File(@"C:\Data\a.raw", "Control"),
                File(@"C:\Data\a.raw", "Treated"),
            };

            // Label-free measures a file once; a repeat is a caller mistake, not a second channel.
            Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.LabelFree(files));
        }

        [Test]
        public void FromSamples_GroupsIsobaricChannelsByFile()
        {
            const string file1 = @"C:\Data\plex1.raw";
            const string file2 = @"C:\Data\plex2.raw";

            var samples = new ISampleInfo[]
            {
                Channel(file1, "126", 126.12776, isReference: true),
                Channel(file1, "127N", 127.12476),
                Channel(file2, "126", 126.12776, isReference: true),
                Channel(file2, "127N", 127.12476),
            };

            var design = SampleExperimentalDesign.FromSamples(samples);

            Assert.Multiple(() =>
            {
                Assert.That(design.FileNameSampleInfoDictionary, Has.Count.EqualTo(2));
                Assert.That(design.FileNameSampleInfoDictionary["plex1.raw"], Has.Length.EqualTo(2));
                Assert.That(design.FileNameSampleInfoDictionary["plex2.raw"], Has.Length.EqualTo(2));
            });
        }

        [Test]
        public void FromSamples_PreservesInputOrderWithinAFile()
        {
            const string path = @"C:\Data\plex.raw";
            var samples = new ISampleInfo[]
            {
                Channel(path, "126", 126.12776),
                Channel(path, "127N", 127.12476),
                Channel(path, "127C", 127.13108),
            };

            var design = SampleExperimentalDesign.FromSamples(samples);

            Assert.That(design.FileNameSampleInfoDictionary["plex.raw"]
                    .Cast<IsobaricQuantSampleInfo>().Select(c => c.ChannelLabel),
                Is.EqualTo(new[] { "126", "127N", "127C" }));
        }

        [Test]
        public void FromSamples_HandlesLabelFreeSamplesToo()
        {
            var samples = new ISampleInfo[]
            {
                File(@"C:\Data\a.raw"),
                File(@"C:\Data\b.raw"),
            };

            var design = SampleExperimentalDesign.FromSamples(samples);

            Assert.Multiple(() =>
            {
                Assert.That(design.FileNameSampleInfoDictionary, Has.Count.EqualTo(2));
                Assert.That(design.FileNameSampleInfoDictionary["a.raw"], Has.Length.EqualTo(1));
            });
        }

        [Test]
        public void FromSamples_RejectsASampleThatNamesNoFile()
        {
            var samples = new ISampleInfo[] { File(string.Empty) };

            var ex = Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.FromSamples(samples));
            Assert.That(ex.Message, Does.Contain("names no file"));
        }

        [Test]
        public void FromSamples_RejectsNull()
        {
            Assert.Throws<ArgumentNullException>(() => SampleExperimentalDesign.FromSamples(null));
        }

        [Test]
        public void LabelFree_RejectsNull()
        {
            Assert.Throws<ArgumentNullException>(() => SampleExperimentalDesign.LabelFree(null));
        }

        [Test]
        public void Add_RejectsANullSampleArray()
        {
            var design = new SampleExperimentalDesign();

            // Distinct from passing no samples at all: params gives an empty array there, null here.
            var ex = Assert.Throws<ArgumentException>(
                () => design.Add("run1.raw", (ISampleInfo[])null));
            Assert.That(ex.Message, Does.Contain("no samples"));
        }

        [Test]
        public void Add_RejectsAPathWithNoFileNameComponent()
        {
            var design = new SampleExperimentalDesign();

            // A directory is not a file, and Path.GetFileName gives back nothing to key on.
            var ex = Assert.Throws<ArgumentException>(
                () => design.Add(@"C:\Data\Experiment\", File(@"C:\Data\Experiment\run1.raw")));
            Assert.That(ex.Message, Does.Contain("no file name component"));
        }

        [Test]
        public void FromSamples_RejectsANullSampleInTheSequence()
        {
            var samples = new ISampleInfo[] { File(@"C:\Data\a.raw"), null };

            var ex = Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.FromSamples(samples));
            Assert.That(ex.Message, Does.Contain("index 1"));
        }

        [Test]
        public void LabelFree_RejectsANullFile()
        {
            var files = new[] { File(@"C:\Data\a.raw"), null };

            var ex = Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.LabelFree(files));
            Assert.That(ex.Message, Does.Contain("null file"));
        }

        [Test]
        public void LabelFree_RejectsAFileWithNoPath()
        {
            var files = new[] { File(string.Empty) };

            var ex = Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.LabelFree(files));
            Assert.That(ex.Message, Does.Contain("no path"));
        }

        [Test]
        public void FromSamples_IsCaseInsensitiveAcrossFiles()
        {
            // Two spellings of one file are one file, matching the lookup the engine will do.
            var samples = new ISampleInfo[]
            {
                Channel(@"C:\Data\Plex.raw", "126", 126.12776),
                Channel(@"C:\Data\plex.raw", "127N", 127.12476),
            };

            var design = SampleExperimentalDesign.FromSamples(samples);

            Assert.That(design.FileNameSampleInfoDictionary, Has.Count.EqualTo(1));
            Assert.That(design.FileNameSampleInfoDictionary["plex.raw"], Has.Length.EqualTo(2));
        }

        [Test]
        public void FromSamples_RejectsTwoDifferentFilesWithTheSameName()
        {
            // Same basename, different directories. Grouping them would put one run's channels in the
            // other run's row, and the design could not tell them apart afterwards.
            var samples = new ISampleInfo[]
            {
                Channel(@"C:\Data\RunA\plex.raw", "126", 126.12776),
                Channel(@"C:\Data\RunB\plex.raw", "127N", 127.12476),
            };

            var ex = Assert.Throws<ArgumentException>(() => SampleExperimentalDesign.FromSamples(samples));
            Assert.Multiple(() =>
            {
                Assert.That(ex.Message, Does.Contain("2 different files"));
                Assert.That(ex.Message, Does.Contain("RunA"));
                Assert.That(ex.Message, Does.Contain("RunB"));
            });
        }

        [Test]
        public void FromSamples_StillGroupsChannelsThatShareOnePath()
        {
            // The case the rejection must not catch: many channels, one file.
            const string path = @"C:\Data\RunA\plex.raw";
            var samples = new ISampleInfo[]
            {
                Channel(path, "126", 126.12776),
                Channel(path, "127N", 127.12476),
                Channel(path, "127C", 127.13108),
            };

            var design = SampleExperimentalDesign.FromSamples(samples);

            Assert.That(design.FileNameSampleInfoDictionary["plex.raw"], Has.Length.EqualTo(3));
        }

        [Test]
        public void ImplementsIExperimentalDesign()
        {
            IExperimentalDesign design = SampleExperimentalDesign.LabelFree(new[] { File(@"C:\Data\a.raw") });

            Assert.That(design.FileNameSampleInfoDictionary, Is.Not.Null);
        }
    }
}
