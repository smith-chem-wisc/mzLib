using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using MassSpectrometry;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// SourceFile.InstrumentSerialNumber and SourceFile.AcquisitionStartTime: which instrument, and
    /// when. Both were in the files all along -- mzML's reader skipped MS:1000529 on purpose and
    /// never read run/@startTimeStamp, and the Thermo reader fetched the instrument data and read
    /// only the model -- so a batch or instrument confound could not be seen from mzLib.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestAcquisitionMetadata
    {
        private static string Data(string relativePath) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, relativePath);

        private static SourceFile Read(string relativePath) =>
            MsDataFileReader.GetDataFile(Data(relativePath)).GetSourceFile();

        [Test]
        [TestCase(@"AveragingTests\TestData\TDYeastFractionMS1.mzML", "EXRFSN20410")]
        [TestCase(@"DatabaseTests\sliced_b6.mzML", "FSN20121")]
        [TestCase(@"FlashLFQ\TestData\sliced-mzml.mzML", "SN03001B")]
        public void MzML_ReadsTheSerialNumberFromTheReferenceableParamGroup(string relativePath, string expected)
        {
            // The ProteoWizard layout: the serial sits beside the model in CommonInstrumentParams.
            Assert.That(Read(relativePath).InstrumentSerialNumber, Is.EqualTo(expected));
        }

        [Test]
        public void MzML_ReadsASerialNumberWrittenInlineOnTheConfiguration()
        {
            Assert.That(Read(@"DataFiles\tiny.pwiz.1.1.mzML").InstrumentSerialNumber, Is.EqualTo("23433"));
        }

        [Test]
        public void MzML_ReportsAPlaceholderSerialVerbatim()
        {
            // It is what the file says. Deciding that "Serial Number N/A" means "unknown" is the
            // consumer's call, not the reader's.
            Assert.That(Read(@"DataFiles\noPrecursorScans.mzML").InstrumentSerialNumber,
                Is.EqualTo("Serial Number N/A"));
        }

        [Test]
        public void MzML_WithNeitherValueReportsNull()
        {
            var sourceFile = Read(@"DataFiles\SmallCalibratibleYeast.mzml");

            Assert.That(sourceFile.InstrumentSerialNumber, Is.Null);
            Assert.That(sourceFile.AcquisitionStartTime, Is.Null);
        }

        [Test]
        public void MzML_StartTimeStampInUtcIsReadAsUtc()
        {
            var time = Read(@"DatabaseTests\sliced_b6.mzML").AcquisitionStartTime;

            Assert.That(time, Is.EqualTo(new DateTime(2018, 12, 13, 20, 0, 48, DateTimeKind.Utc)));
            Assert.That(time!.Value.Kind, Is.EqualTo(DateTimeKind.Utc));
        }

        [Test]
        public void MzML_StartTimeStampKeepsFractionalSeconds()
        {
            var time = Read(@"DataFiles\badScan7192.mzML").AcquisitionStartTime;

            Assert.That(time, Is.EqualTo(new DateTime(2021, 1, 13, 9, 34, 11, 285, DateTimeKind.Utc)));
        }

        [Test]
        public void MzML_StartTimeStampWithNoOffsetStaysUnspecified()
        {
            // No offset, no instant: the reader must not guess one from the machine it runs on.
            var time = Read(@"DataFiles\tiny.pwiz.1.1.mzML").AcquisitionStartTime;

            Assert.That(time!.Value.Kind, Is.EqualTo(DateTimeKind.Unspecified));
            Assert.That(time.Value, Is.EqualTo(new DateTime(2007, 6, 27, 15, 23, 45).AddTicks(3500)));
        }

        [Test]
        public void MzML_StartTimeStampWithAnOffsetIsConvertedToUtc_WhateverTheMachineTimeZone()
        {
            // XmlSerializer turns an explicit offset into the READING machine's Local time. Converted
            // back, the instant is the same on every machine.
            string text = File.ReadAllText(Data(@"DataFiles\tiny.pwiz.1.1.mzML"))
                .Replace("startTimeStamp=\"2007-06-27T15:23:45.00035\"", "startTimeStamp=\"2007-06-27T15:23:45+02:00\"");
            string path = Data($"offset_{Guid.NewGuid():N}.mzML");
            File.WriteAllText(path, text);
            try
            {
                var time = MsDataFileReader.GetDataFile(path).GetSourceFile().AcquisitionStartTime;

                Assert.That(time, Is.EqualTo(new DateTime(2007, 6, 27, 13, 23, 45, DateTimeKind.Utc)));
                Assert.That(time!.Value.Kind, Is.EqualTo(DateTimeKind.Utc));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        [TestCase(@"DataFiles\sliced_ethcd.raw", @"DataFiles\sliced_ethcd.mzML", "FSN10189")]
        [TestCase(@"FlashLFQ\TestData\sliced-raw.raw", @"FlashLFQ\TestData\sliced-mzml.mzML", "SN03001B")]
        public void ThermoRaw_SerialNumberMatchesProteoWizardsConversionOfTheSameFile(
            string raw, string mzml, string expected)
        {
            // Each mzML records the SHA-1 of this exact RAW, so it is ProteoWizard's reading of it.
            Assert.That(Read(raw).InstrumentSerialNumber, Is.EqualTo(expected));
            Assert.That(Read(mzml).InstrumentSerialNumber, Is.EqualTo(expected));
        }

        [Test]
        [TestCase(@"DataFiles\sliced_ethcd.raw", 2021, 3, 16, 12, 9, 7)]
        [TestCase(@"FlashLFQ\TestData\sliced-raw.raw", 2017, 7, 11, 9, 55, 53)]
        public void ThermoRaw_StartTimeIsTheAcquisitionComputersWallClock_WithNoKind(
            string raw, int year, int month, int day, int hour, int minute, int second)
        {
            // The Thermo library labels this Utc, but it is local to the acquisition computer:
            // ProteoWizard's mzML of each file says five hours LATER, with "Z" (US Central daylight
            // time, converted on a machine in that zone). The RAW records no offset, so neither
            // instant can be confirmed, and the Kind is cleared.
            var time = Read(raw).AcquisitionStartTime;

            Assert.That(time!.Value.Kind, Is.EqualTo(DateTimeKind.Unspecified));
            Assert.That(new DateTime(time.Value.Year, time.Value.Month, time.Value.Day,
                    time.Value.Hour, time.Value.Minute, time.Value.Second),
                Is.EqualTo(new DateTime(year, month, day, hour, minute, second)));
        }

        [Test]
        public void ThermoRaw_ThatCannotBeOpenedReportsNeitherValue()
        {
            string path = Data($"not_really_{Guid.NewGuid():N}.raw");
            File.WriteAllText(path, "this is not a RAW file");
            try
            {
                var sourceFile = MsDataFileReader.GetDataFile(path).GetSourceFile();

                Assert.That(sourceFile.CheckSum, Is.Not.Empty);
                Assert.That(sourceFile.InstrumentSerialNumber, Is.Null);
                Assert.That(sourceFile.AcquisitionStartTime, Is.Null);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        [TestCase(null, null)]
        [TestCase("", null)]
        [TestCase("   ", null)]
        [TestCase(" FSN10189 ", "FSN10189")]
        public void ThermoRaw_SerialNumberIsTrimmed_AndBlankIsNull(string given, string expected)
        {
            Assert.That(ThermoRawFileReader.BuildSerialNumber(given), Is.EqualTo(expected));
        }

        [Test]
        public void ThermoRaw_StartTimeClearsTheKind_AndAnUnsetDateIsNull()
        {
            var labelledUtc = new DateTime(2021, 3, 16, 12, 9, 7, DateTimeKind.Utc);

            var time = ThermoRawFileReader.BuildAcquisitionStartTime(labelledUtc);

            Assert.That(time!.Value.Kind, Is.EqualTo(DateTimeKind.Unspecified));
            Assert.That(time.Value.Ticks, Is.EqualTo(labelledUtc.Ticks), "the wall-clock value is kept");
            Assert.That(ThermoRawFileReader.BuildAcquisitionStartTime(null), Is.Null);
            Assert.That(ThermoRawFileReader.BuildAcquisitionStartTime(default(DateTime)), Is.Null);
        }

        [Test]
        [TestCase(@"DatabaseTests\sliced_b6.mzML")]
        [TestCase(@"DataFiles\sliced_ethcd.raw")]
        public void BothValues_SurviveAWrite(string relativePath)
        {
            // Calibration and averaging rewrite the mzML; the serial and the time, with its Kind,
            // must come back out as they went in, whichever reader produced them. The mzML case
            // carries a Utc time and the RAW case an Unspecified one.
            var original = MsDataFileReader.GetDataFile(Data(relativePath));
            original.LoadAllStaticData();
            var before = original.GetSourceFile();
            Assert.That(before.InstrumentSerialNumber, Is.Not.Null, "precondition");
            Assert.That(before.AcquisitionStartTime, Is.Not.Null, "precondition");

            string outputPath = Data($"acquisition_rt_{Guid.NewGuid():N}.mzML");
            try
            {
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(original, outputPath, false);
                var after = MsDataFileReader.GetDataFile(outputPath).GetSourceFile();

                Assert.That(after.InstrumentSerialNumber, Is.EqualTo(before.InstrumentSerialNumber));
                Assert.That(after.AcquisitionStartTime, Is.EqualTo(before.AcquisitionStartTime));
                Assert.That(after.AcquisitionStartTime!.Value.Kind, Is.EqualTo(before.AcquisitionStartTime!.Value.Kind));
                Assert.That(after.InstrumentModel?.Accession, Is.Not.EqualTo("MS:1000529"),
                    "the serial written beside the model must not be read back as the model");
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void UnknownValues_AreNotWritten()
        {
            var original = MsDataFileReader.GetDataFile(Data(@"DataFiles\SmallCalibratibleYeast.mzml"));
            original.LoadAllStaticData();

            string outputPath = Data($"acquisition_none_{Guid.NewGuid():N}.mzML");
            try
            {
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(original, outputPath, false);
                string text = File.ReadAllText(outputPath);

                Assert.That(text, Does.Not.Contain("MS:1000529"));
                Assert.That(text, Does.Not.Contain("startTimeStamp"));
            }
            finally
            {
                if (File.Exists(outputPath)) File.Delete(outputPath);
            }
        }

        [Test]
        public void BothValues_SurviveASnip()
        {
            // ExportSnipAsMzML copies SourceFile field by field, so a new init-only property is
            // silently dropped there unless it is carried explicitly.
            var original = MsDataFileReader.GetDataFile(Data(@"DatabaseTests\sliced_b6.mzML"));
            original.LoadAllStaticData();
            var before = original.GetSourceFile();

            string snipPath = original.ExportSnipAsMzML(1, 5);
            try
            {
                var after = MsDataFileReader.GetDataFile(snipPath).GetSourceFile();

                Assert.That(after.InstrumentSerialNumber, Is.EqualTo(before.InstrumentSerialNumber));
                Assert.That(after.AcquisitionStartTime, Is.EqualTo(before.AcquisitionStartTime));
            }
            finally
            {
                if (File.Exists(snipPath)) File.Delete(snipPath);
            }
        }

        [Test]
        public void BothValues_DefaultToNullWhenUnset()
        {
            var sourceFile = new SourceFile("no nativeID format", "mzML format", "abc", "SHA-1", "id");

            Assert.That(sourceFile.InstrumentSerialNumber, Is.Null);
            Assert.That(sourceFile.AcquisitionStartTime, Is.Null);
        }
    }
}
