using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using SpectralAveraging;

namespace Test.AveragingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestAveragingSpectraWriteFile
    {
        public static string OutputDirectory;
        public static string SpectraPath;
        public static SpectralAveragingParameters Parameters;
        public static List<MsDataScan> Scans;
        public static MsDataScan[] DdaCompositeSpectra;

        [OneTimeSetUp]
        public static void OneTimeSetup()
        {
            Parameters = new SpectralAveragingParameters();
            OutputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"AveragingTestData");
            SpectraPath = Path.Combine(OutputDirectory, "TDYeastFractionMS1.mzML");
            Scans = MsDataFileReader.GetDataFile(SpectraPath).GetAllScansList().Take(50).ToList();

            Parameters.SpectraFileAveragingType = SpectraFileAveragingType.AverageDdaScansWithOverlap;
            DdaCompositeSpectra = SpectraFileAveraging.AverageSpectraFile(Scans, Parameters);
            Assert.That(DdaCompositeSpectra.Length > 1);
        }

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            if (Directory.Exists(Path.Combine(OutputDirectory, "AveragedSpectra")))
                Directory.Delete(Path.Combine(OutputDirectory, "AveragedSpectra"), true);
            string[] filesToDelete =
                Directory.GetFiles(OutputDirectory).Where(p => p.Contains("Averaged", StringComparison.InvariantCultureIgnoreCase)).ToArray();
            foreach (var file in filesToDelete)
            {
                File.Delete(file);
            }
            Assert.That(Directory.GetFiles(OutputDirectory).Length == 2);

            foreach (var directoryToDelete in Directory.GetDirectories(OutputDirectory).Where(p => p.Contains("Averaged", StringComparison.InvariantCultureIgnoreCase)))
            {
                Directory.Delete(directoryToDelete, true);
            }
            Assert.That(!Directory.GetDirectories(OutputDirectory).Any());
        }

        [Test]
        public static void OutputAveragedSpectraAsMzMLTest()
        {
            // test that it outputs correctly
            Assert.That(Parameters.OutputType == OutputType.MzML);
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath);
            string averagedSpectraPath = Path.Combine(OutputDirectory,
                Path.GetFileNameWithoutExtension(SpectraPath) + "-averaged.mzML");
            Assert.That(File.Exists(averagedSpectraPath));

            MsDataScan[] loadedScans = MsDataFileReader.GetDataFile(averagedSpectraPath).GetAllScansList().ToArray();
            for (var i = 0; i < loadedScans.Length; i++)
            {
                for (int j = 0; j < loadedScans[i].MassSpectrum.YArray.Length; j++)
                {
                    Assert.That(Math.Abs(loadedScans[i].MassSpectrum.XArray[j] - DdaCompositeSpectra[i].MassSpectrum.XArray[j]) < 0.0001);
                    Assert.That(Math.Abs(loadedScans[i].MassSpectrum.YArray[j] - DdaCompositeSpectra[i].MassSpectrum.YArray[j]) < 0.0001);
                }
            }

            // test errors
            var exception = Assert.Throws<MzLibException>(() =>
            {
                AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, "");
            });
            Assert.That(exception.Message == "Cannot Access Spectra Directory");
        }

        [Test]
        public static void TestRepeatOutputToSameDirectoryWithSameNameMzML()
        {
            string averagedSpectraPath = Path.Combine(OutputDirectory, Path.GetFileNameWithoutExtension(SpectraPath) + "-averaged.mzML");
            string secondOutputPath = Path.Combine(OutputDirectory, Path.GetFileNameWithoutExtension(SpectraPath) + "-averaged(1).mzML");
            string thirdOutputPath = Path.Combine(OutputDirectory, Path.GetFileNameWithoutExtension(SpectraPath) + "-averaged(2).mzML");

            // clean up from any previous tests
            foreach (var averagedFile in Directory.GetFiles(OutputDirectory)
                         .Where(p => p.Contains("Averaged", StringComparison.InvariantCultureIgnoreCase)))
            {
                File.Delete(averagedFile);
            }

            // write three times and ensure unique outputs are generated
            Parameters.OutputType = OutputType.MzML;
            Assert.That(!File.Exists(averagedSpectraPath));
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath);
            Assert.That(File.Exists(averagedSpectraPath));
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath);
            Assert.That(File.Exists(secondOutputPath));
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath);
            Assert.That(File.Exists(thirdOutputPath));
        }

        [Test]
        public static void TestOutputToCustomDirectoryAndNameMzML()
        {
            // output to a different directory than the files were originally in
            Parameters.OutputType = OutputType.MzML;
            string customDestinationDirectory = Path.Combine(OutputDirectory, "NewTestingDirectory");
            string customDestinationDirectory2 = Path.Combine(OutputDirectory, "NewTestingDirectory2");
            Directory.CreateDirectory(customDestinationDirectory);
            string customName = "AveragedSpectra";

            // custom destination, original name
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath,
                customDestinationDirectory);
            var files = Directory.GetFiles(customDestinationDirectory);
            Assert.That(files.Length == 1);
            Assert.That(files.Contains(Path.Combine(customDestinationDirectory,
                Path.GetFileNameWithoutExtension(SpectraPath) + "-averaged.mzML")));

            // custom destination, custom name
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath,
                customDestinationDirectory, customName);
            files = Directory.GetFiles(customDestinationDirectory);
            Assert.That(files.Length == 2);
            Assert.That(files.Contains(Path.Combine(customDestinationDirectory, customName + ".mzML")));
            
            // custom destination, custom name : directory not created before run
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath,
                customDestinationDirectory2, customName);
            files = Directory.GetFiles(customDestinationDirectory2);
            Assert.That(files.Length == 1);
            Assert.That(files.Contains(Path.Combine(customDestinationDirectory2, customName + ".mzML")));

            // original destination, custom name
            AveragedSpectraWriter.WriteAveragedScans(DdaCompositeSpectra, Parameters, SpectraPath, null, customName);
            Assert.That(File.Exists(Path.Combine(OutputDirectory, customName + ".mzML")));

            // clean up custom destinations
            Directory.Delete(customDestinationDirectory, true);
            Directory.Delete(customDestinationDirectory2, true);
        }

       

    }
}
