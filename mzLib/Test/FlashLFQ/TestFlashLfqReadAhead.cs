using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.FlashLFQ
{
    /// <summary>
    /// FlashLFQ reads and indexes the next spectra file while the current one is quantified, when there is a thread
    /// to spare and memory for both indexes. The index of a file depends only on that file, so the results must be
    /// exactly those of reading each file in turn, and a file that fails to read must fail where it did before.
    /// </summary>
    [TestFixture]
    public static class TestFlashLfqReadAhead
    {
        private static string DataDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData");

        [Test]
        public static void ReadingAheadGivesTheSameResults()
        {
            string work = Path.Combine(TestContext.CurrentContext.WorkDirectory, nameof(ReadingAheadGivesTheSameResults));
            if (Directory.Exists(work)) Directory.Delete(work, true);
            Directory.CreateDirectory(work);

            // Three files, so one read-ahead follows another: the two K562 files and a copy of the first under another name
            string[] names = { "20100614_Velos1_TaGe_SA_K562_3", "20100614_Velos1_TaGe_SA_K562_4", "20100614_Velos1_TaGe_SA_K562_3_copy" };
            var files = new List<SpectraFileInfo>();
            for (int i = 0; i < names.Length; i++)
            {
                string source = Path.Combine(DataDirectory, names[i].Replace("_copy", "") + ".mzML");
                string path = Path.Combine(work, names[i] + ".mzML");
                File.Copy(source, path, true);
                files.Add(new SpectraFileInfo(path, "a", i, 0, 0));
            }

            GC.Collect(); // the memory check needs one collection's worth of numbers
            string sequential = Quantify(files, maxThreads: 1, out int readAheadSequential);
            string readingAhead = Quantify(files, maxThreads: 4, out int readAheadParallel);

            Assert.That(readAheadSequential, Is.EqualTo(0), "a single-threaded run must not read ahead");
            Assert.That(readingAhead, Is.EqualTo(sequential));
            Assume.That(readAheadParallel, Is.EqualTo(files.Count - 1), "this machine had no memory to spare for reading ahead");

            Directory.Delete(work, true);
        }

        [Test]
        public static void AFileThatCannotBeReadFailsAsBefore()
        {
            string work = Path.Combine(TestContext.CurrentContext.WorkDirectory, nameof(AFileThatCannotBeReadFailsAsBefore));
            if (Directory.Exists(work)) Directory.Delete(work, true);
            Directory.CreateDirectory(work);

            string good = Path.Combine(work, "20100614_Velos1_TaGe_SA_K562_3.mzML");
            File.Copy(Path.Combine(DataDirectory, "20100614_Velos1_TaGe_SA_K562_3.mzML"), good, true);
            string bad = Path.Combine(work, "20100614_Velos1_TaGe_SA_K562_4.mzML");
            File.WriteAllText(bad, string.Concat(Enumerable.Repeat("this is not an mzML file\n", 100_000)));
            var files = new List<SpectraFileInfo> { new SpectraFileInfo(good, "a", 0, 0, 0), new SpectraFileInfo(bad, "a", 1, 0, 0) };

            GC.Collect();
            string sequential = Outcome(files, maxThreads: 1);
            string readingAhead = Outcome(files, maxThreads: 4);
            Assert.That(readingAhead, Is.EqualTo(sequential));

            Directory.Delete(work, true);
        }

        private static string Outcome(List<SpectraFileInfo> files, int maxThreads)
        {
            try
            {
                return Quantify(files, maxThreads, out _);
            }
            catch (Exception e)
            {
                return "threw " + e.GetType().FullName + ": " + e.Message;
            }
        }

        private static string Quantify(List<SpectraFileInfo> files, int maxThreads, out int filesReadAhead)
        {
            var engine = new FlashLfqEngine(ReadIdentifications(files), matchBetweenRuns: true, normalize: true, maxThreads: maxThreads,
                silent: true, randomSeed: 42, matchBetweenRunsFdrThreshold: 0.5);
            var results = engine.Run();
            filesReadAhead = engine.FilesReadAhead;

            string folder = Path.GetDirectoryName(files[0].FullFilePathWithExtension);
            string peaks = Path.Combine(folder, $"peaks_{maxThreads}.tsv");
            string peptides = Path.Combine(folder, $"peptides_{maxThreads}.tsv");
            string proteins = Path.Combine(folder, $"proteins_{maxThreads}.tsv");
            results.WriteResults(peaks, peptides, proteins, null, true);
            return File.ReadAllText(peaks) + File.ReadAllText(peptides) + File.ReadAllText(proteins);
        }

        private static List<Identification> ReadIdentifications(List<SpectraFileInfo> files)
        {
            var identifications = new List<Identification>();
            foreach (var psm in SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(DataDirectory, "AllPSMs.psmtsv"), out _))
            {
                foreach (SpectraFileInfo file in files.Where(f => f.FilenameWithoutExtension.Replace("_copy", "") == psm.FileNameWithoutExtension))
                {
                    var proteins = psm.ProteinAccession.Split('|').Select(a => new ProteinGroup(a, "", "")).ToList();
                    identifications.Add(new Identification(file, psm.BaseSeq, psm.FullSequence,
                        (double)psm.MonoisotopicMass, (double)psm.RetentionTime, psm.PrecursorCharge, proteins,
                        decoy: psm.DecoyContamTarget == "D"));
                }
            }
            return identifications;
        }
    }
}
