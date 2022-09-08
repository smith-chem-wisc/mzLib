using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using IO.MzML;
using MassSpectrometry;
using NUnit.Framework;
using SpectralAveragingExtensions;
using OutputType = SpectralAveragingExtensions.OutputType;
using SpectraFileProcessingType = SpectralAveragingExtensions.SpectraFileProcessingType;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class AveragingSpectraOutputFileTests
    {
        public static string OutputDirectory;
        public static string SpectraPath;
        public static MzLibSpectralAveragingOptions Options;
        public static List<MsDataScan> Scans;
        public static MsDataScan[] DdaCompositeSpectra;
        public static MsDataScan[] AverageAllCompositeSpectra;

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Options = new MzLibSpectralAveragingOptions();
            OutputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"AveragingTestData");
            SpectraPath = Path.Combine(OutputDirectory, "TDYeastFractionMS1.mzML");
            Scans = SpectraFileHandler.LoadAllScansFromFile(SpectraPath);
            Options.SpectraFileProcessingType = SpectraFileProcessingType.AverageDDAScansWithOverlap;
            DdaCompositeSpectra = SpectraFileProcessing.ProcessSpectra(Scans, Options);
            Assert.That(DdaCompositeSpectra.Length > 1);
            Options.SpectraFileProcessingType = SpectraFileProcessingType.AverageAll;
            AverageAllCompositeSpectra = SpectraFileProcessing.ProcessSpectra(Scans, Options);
        }

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            if (Directory.Exists(Path.Combine(OutputDirectory, "AveragedSpectra")))
                Directory.Delete(Path.Combine(OutputDirectory, "AveragedSpectra"), true);
            string[] filesToDelete =
                Directory.GetFiles(OutputDirectory).Where(p => p.Contains("Averaged")).ToArray();
            foreach (var file in filesToDelete)
            {
                File.Delete(file);
            }
            Assert.That(Directory.GetFiles(OutputDirectory).Length == 2);
        }

        [Test]
        public static void OutputAveragedSpectraAsMzMLTest()
        {
            Assert.That(Options.OutputType == OutputType.mzML);
            AveragedSpectraOutputter.OutputAveragedScans(DdaCompositeSpectra, Options, SpectraPath);
            string averagedSpectraPath = Path.Combine(OutputDirectory,
                "Averaged_" + Path.GetFileNameWithoutExtension(SpectraPath) + ".mzML");
            Assert.That(File.Exists(averagedSpectraPath));

            MsDataScan[] loadedScans = Mzml.LoadAllStaticData(averagedSpectraPath).GetAllScansList().ToArray();
            for (var i = 0; i < loadedScans.Length; i++)
            {
                Dictionary<double, double> mzandInt = DdaCompositeSpectra[i].MassSpectrum.XArray
                    .Zip(DdaCompositeSpectra[i].MassSpectrum.YArray, (x, y) => new { x, y })
                    .ToDictionary(p => p.x, p => p.y);
                Dictionary<double, double> trimmedMzAndInt =
                    mzandInt.Where(p => p.Value != 0).ToDictionary(p => p.Key, p => p.Value);
                MzSpectrum trimmedDdaSpectrum = new MzSpectrum(trimmedMzAndInt.Keys.ToArray(), trimmedMzAndInt.Values.ToArray(), false);
                Assert.That(loadedScans[i].MassSpectrum.Equals(trimmedDdaSpectrum));
            }
        }

        [Test]
        public static void OutputAveragedSpectraAsTxt()
        {
            Options.OutputType = OutputType.txt;
            Assert.That(Options.OutputType == OutputType.txt);

            // test outputs where there will be multiple spectra
            Options.SpectraFileProcessingType = SpectraFileProcessingType.AverageDDAScansWithOverlap;
            AveragedSpectraOutputter.OutputAveragedScans(DdaCompositeSpectra, Options, SpectraPath);
            Assert.That(Directory.Exists(Path.Combine(OutputDirectory, "AveragedSpectra")));
            string[] txtFiles = Directory.GetFiles(Path.Combine(OutputDirectory, "AveragedSpectra"))
                .OrderBy(p => double.Parse(p.Split('_').Last().Replace(".txt", ""))).ToArray();
            Assert.That(txtFiles.Length == DdaCompositeSpectra.Length);
            for (var i = 0; i < txtFiles.Length; i++)
            {

                double[] xArray = new double[DdaCompositeSpectra[i].MassSpectrum.XArray.Length];
                double[] yArray = new double[DdaCompositeSpectra[i].MassSpectrum.YArray.Length];
                string[] mzAndInts = File.ReadAllLines(txtFiles[i]);
                for (var j = 0; j < mzAndInts.Length; j++)
                {
                    xArray[j] = double.Parse(mzAndInts[j].Split(',')[0]);
                    yArray[j] = double.Parse(mzAndInts[j].Split(',')[1]);
                }

                MzSpectrum loadedSpectra = new(xArray, yArray, false);
                Assert.That(loadedSpectra.Equals(DdaCompositeSpectra[i].MassSpectrum));
            }

            // test average all output
            Options.SpectraFileProcessingType = SpectraFileProcessingType.AverageAll;
            AveragedSpectraOutputter.OutputAveragedScans(AverageAllCompositeSpectra, Options, SpectraPath);
            string averagedSpectraPath = Path.Combine(OutputDirectory,
                "Averaged_" + Path.GetFileNameWithoutExtension(SpectraPath) + ".txt");
            Assert.That(File.Exists(averagedSpectraPath));
            string[] mzAndInt = File.ReadAllLines(averagedSpectraPath);
            double[] xArr = mzAndInt.Select(p => double.Parse(p.Split(',')[0])).ToArray();
            double[] yArr = mzAndInt.Select(p => double.Parse(p.Split(',')[1])).ToArray();
            MzSpectrum loadedSpectrum = new(xArr, yArr, true);
            Assert.That(loadedSpectrum.Equals(AverageAllCompositeSpectra[0].MassSpectrum));
        }
    }
}
