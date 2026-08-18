using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using Generated = Readers.Generated;

namespace Test.FileReadingTests
{
    /// <summary>
    /// The rules that decide WHICH cvParam is the instrument model, driven over inputs no committed
    /// file provides: a run that names a default configuration, a param group crowded with the
    /// terms that sit beside a model in real files, and a configuration that states the model
    /// itself while also naming a group.
    ///
    /// Every real converted file states the model the same way, so on committed fixtures these
    /// rules are only ever driven down their accepting side -- a reader that returned the first
    /// MS: term it saw would pass all of them, and would report a serial number as the instrument.
    ///
    /// The input is built as a Generated.mzMLType object graph and serialized with the writer's own
    /// serializer rather than written as XML text, so a schema change cannot leave a hand-written
    /// fixture silently unparsed. Same shape as the mzIdentML round trips in TestMzML.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestInstrumentModelSelection
    {
        private const string LumosAccession = "MS:1002732";
        private const string LumosName = "Orbitrap Fusion Lumos";
        private const string VelosAccession = "MS:1001742";
        private const string VelosName = "LTQ Orbitrap Velos";

        /// <summary>
        /// Serializes the graph to a temporary mzML, reads it back through the production reader,
        /// and returns the instrument model it found.
        /// </summary>
        private static CvParam ReadInstrumentModel(Generated.mzMLType mzML)
        {
            string path = Path.Combine(
                TestContext.CurrentContext.TestDirectory, $"instrument_edge_{Guid.NewGuid():N}.mzML");
            try
            {
                using (var writer = new StreamWriter(path))
                {
                    MzmlMethods.mzmlSerializer.Serialize(writer, mzML);
                }

                // GetSourceFile is the whole route to the reader's instrument-model logic; the scan
                // list is deliberately absent, so LoadAllStaticData is neither needed nor valid here.
                return new Mzml(path).GetSourceFile().InstrumentModel;
            }
            finally
            {
                if (File.Exists(path))
                {
                    File.Delete(path);
                }
            }
        }

        private static Generated.CVParamType Cv(string accession, string name, string value = "") =>
            new Generated.CVParamType { cvRef = "MS", accession = accession, name = name, value = value };

        private static Generated.ReferenceableParamGroupType Group(
            string id, params Generated.CVParamType[] cvParams) =>
            new Generated.ReferenceableParamGroupType { id = id, cvParam = cvParams };

        /// <summary>
        /// An mzML carrying nothing but the instrument metadata under test. fileDescription is
        /// required -- GetSourceFile reads it before anything else -- and an absent sourceFileList
        /// sends it down its checksum branch, which is what a file written by hand would do.
        /// </summary>
        private static Generated.mzMLType FileWith(
            Generated.InstrumentConfigurationType[] configurations,
            Generated.ReferenceableParamGroupType[] groups = null,
            string defaultConfigurationRef = null)
        {
            var mzML = new Generated.mzMLType
            {
                fileDescription = new Generated.FileDescriptionType()
            };

            if (configurations != null)
            {
                mzML.instrumentConfigurationList = new Generated.InstrumentConfigurationListType
                {
                    count = configurations.Length.ToString(),
                    instrumentConfiguration = configurations
                };
            }

            if (groups != null)
            {
                mzML.referenceableParamGroupList = new Generated.ReferenceableParamGroupListType
                {
                    count = groups.Length.ToString(),
                    referenceableParamGroup = groups
                };
            }

            if (defaultConfigurationRef != null)
            {
                mzML.run = new Generated.RunType
                {
                    id = "run1",
                    defaultInstrumentConfigurationRef = defaultConfigurationRef
                };
            }

            return mzML;
        }

        [Test]
        public void MzML_ReadsTheInstrumentTheRunDeclaresAsDefault()
        {
            // Instrument model is a property of the run, not of a scan. A file that lists several
            // configurations and names one as the run's default is stating which instrument
            // produced the data; taking the first entry instead would report the wrong instrument
            // on any file whose default is not listed first, and report it confidently.
            var model = ReadInstrumentModel(FileWith(
                new[]
                {
                    new Generated.InstrumentConfigurationType
                    { id = "IC1", cvParam = new[] { Cv(VelosAccession, VelosName) } },
                    new Generated.InstrumentConfigurationType
                    { id = "IC2", cvParam = new[] { Cv(LumosAccession, LumosName) } }
                },
                defaultConfigurationRef: "IC2"));

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Accession, Is.EqualTo(LumosAccession));
            Assert.That(model.Name, Is.EqualTo(LumosName));
        }

        [Test]
        public void MzML_FallsBackToTheFirstConfigurationWhenTheRunNamesNoDefault()
        {
            // defaultInstrumentConfigurationRef is required by the schema but absent in practice;
            // the first configuration is then the only defensible answer, and is still an answer.
            var model = ReadInstrumentModel(FileWith(
                new[]
                {
                    new Generated.InstrumentConfigurationType
                    { id = "IC1", cvParam = new[] { Cv(LumosAccession, LumosName) } },
                    new Generated.InstrumentConfigurationType
                    { id = "IC2", cvParam = new[] { Cv(VelosAccession, VelosName) } }
                }));

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Accession, Is.EqualTo(LumosAccession));
        }

        [Test]
        public void MzML_PicksTheModelOutOfAGroupCrowdedWithTermsThatAreNotModels()
        {
            // The whole reason the model cannot be identified by a list of accessions: it is
            // surrounded by terms that look just like it. Every decoy here is one the rules reject
            // for a different reason, and the model is deliberately last, so a reader that took the
            // first MS: term, or the first valueless one, or the first named one, fails.
            //
            // The file also declares the spectrum-parameter group pwiz always writes, listed first,
            // so resolving the configuration's reference has to walk past a group it does not name.
            var model = ReadInstrumentModel(FileWith(
                new[]
                {
                    new Generated.InstrumentConfigurationType
                    {
                        id = "IC1",
                        referenceableParamGroupRef = new[]
                        {
                            new Generated.ReferenceableParamGroupRefType { @ref = "CommonInstrumentParams" }
                        }
                    }
                },
                new[]
                {
                    Group("SpectrumParams",
                        Cv("MS:1000579", "MS1 spectrum"),
                        Cv("MS:1000130", "positive scan")),
                    Group("CommonInstrumentParams",
                        Cv("UO:0000010", "second"),                       // another vocabulary
                        Cv(null, "term with no accession"),               // unresolvable
                        Cv("MS:1000529", "instrument serial number", "EXRFSN20410"), // beside it in every pwiz file
                        Cv("MS:1009999", "an attribute this reader has never seen", "42"), // unlisted, but valued
                        Cv("MS:1009998", "   "),                          // accession with no name
                        Cv(LumosAccession, LumosName))
                }));

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Accession, Is.EqualTo(LumosAccession));
            Assert.That(model.Name, Is.EqualTo(LumosName));
            Assert.That(model.Value, Is.Empty, "an instrument model is a bare presence flag");
        }

        [Test]
        public void MzML_InstrumentStatedOnTheConfigurationWinsOverTheSharedGroup()
        {
            // A group is shared across configurations; a term on the configuration itself is
            // specific to it. When the two disagree the specific one is the answer -- otherwise a
            // multi-instrument file reports the same model for every configuration in it.
            var model = ReadInstrumentModel(FileWith(
                new[]
                {
                    new Generated.InstrumentConfigurationType
                    {
                        id = "IC1",
                        cvParam = new[] { Cv(LumosAccession, LumosName) },
                        referenceableParamGroupRef = new[]
                        {
                            new Generated.ReferenceableParamGroupRefType { @ref = "CommonInstrumentParams" }
                        }
                    }
                },
                new[] { Group("CommonInstrumentParams", Cv(VelosAccession, VelosName)) }));

            Assert.That(model, Is.Not.Null);
            Assert.That(model.Name, Is.EqualTo(LumosName));
        }

        [Test]
        public void MzML_ThatDeclaresNoInstrumentAtAllReportsNoModel()
        {
            // Reporting nothing is the correct answer, and it must be reached without throwing:
            // the instrument is metadata, and a file is still readable without it.
            Assert.That(ReadInstrumentModel(FileWith(null)), Is.Null,
                "no instrumentConfigurationList");

            Assert.That(
                ReadInstrumentModel(FileWith(Array.Empty<Generated.InstrumentConfigurationType>())),
                Is.Null, "an empty instrumentConfigurationList");
        }
    }
}
