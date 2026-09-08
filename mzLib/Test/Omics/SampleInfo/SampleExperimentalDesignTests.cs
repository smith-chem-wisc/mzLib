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
