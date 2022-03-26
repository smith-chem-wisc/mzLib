using MzIdentML;
using NUnit.Framework;
using System;
using System.IO;
using System.Text;
using System.Xml;
using System.Xml.Serialization;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestMzIdentML
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

        /// <summary>
        /// Basic test covered by codecov
        /// </summary>
        [Test]
        public static void WriteMzID110Test()
        {
            string version = "1.1.0";
            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType110()
            {
                version = version,
                id = "",
            };

            _mzid.Provider = new mzIdentML110.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML110.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML110.Generated.RoleType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            XmlWriter writer = XmlWriter.Create(Path.Combine(Environment.CurrentDirectory, $"{version}filename"), settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        /// <summary>
        /// Basic test covered by codecov
        /// </summary>
        [Test]
        public static void WriteMzID111Test()
        {
            string version = "1.1.1";
            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML111.Generated.MzIdentMLType111));
            var _mzid = new mzIdentML111.Generated.MzIdentMLType111()
            {
                version = version,
                id = "",
            };

            _mzid.Provider = new mzIdentML111.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML111.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML111.Generated.RoleType()
                    {
                        cvParam = new mzIdentML111.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            XmlWriter writer = XmlWriter.Create(Path.Combine(Environment.CurrentDirectory, $"{version}filename"), settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        /// <summary>
        /// Basic test covered by codecov
        /// </summary>
        [Test]
        public static void WriteMzID120Test()
        {
            string version = "1.2.0";
            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
            var _mzid = new mzIdentML120.Generated.MzIdentMLType120()
            {
                version = version,
                id = "",
            };

            _mzid.Provider = new mzIdentML120.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML120.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML120.Generated.RoleType()
                    {
                        cvParam = new mzIdentML120.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            XmlWriter writer = XmlWriter.Create(Path.Combine(Environment.CurrentDirectory, $"{version}filename"), settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();
        }

        [Test]
        public static void Mzid111Test()
        {
            string version = "1.1.1";
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML111.Generated.MzIdentMLType111));
            var _mzid = new mzIdentML111.Generated.MzIdentMLType111
            {
                version = version,
                id = "",
                DataCollection = new mzIdentML111.Generated.DataCollectionType()
            };
            _mzid.Provider = new mzIdentML111.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML111.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML111.Generated.RoleType()
                    {
                        cvParam = new mzIdentML111.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML111.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML111.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML111.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML111.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML111.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML111.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML111.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML111.Generated.CVParamType[]
                {
                    new mzIdentML111.Generated.CVParamType
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] =
                new mzIdentML111.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation =
                new mzIdentML111.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] =
                new mzIdentML111.Generated.IonTypeType
                {
                    FragmentArray = new mzIdentML111.Generated.FragmentArrayType[1]
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] =
                new mzIdentML111.Generated.FragmentArrayType
                {
                    values = new float[3] { 200, 300, 400 }
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef =
                new mzIdentML111.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] =
                new mzIdentML111.Generated.PeptideEvidenceRefType
                {
                    peptideEvidence_ref = "PE_1"
                };
            _mzid.DataCollection.Inputs = new mzIdentML111.Generated.InputsType
            {
                SpectraData = new mzIdentML111.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML111.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML111.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML111.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML111.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML111.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML111.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML111.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML111.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML111.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML111.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML111.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML111.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML111.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML111.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML111.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML111.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML111.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };

            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlWriter writer = XmlWriter.Create("myIdentifications.mzid", settings);

            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public static void Mzid120Test()
        {
            string version = "1.2.0";
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
            var _mzid = new mzIdentML120.Generated.MzIdentMLType120
            {
                version = version,
                id = "",
                DataCollection = new mzIdentML120.Generated.DataCollectionType()
            };
            _mzid.Provider = new mzIdentML120.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML120.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML120.Generated.RoleType()
                    {
                        cvParam = new mzIdentML120.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            _mzid.DataCollection.AnalysisData = new mzIdentML120.Generated.AnalysisDataType()
            {
                SpectrumIdentificationList = new mzIdentML120.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML120.Generated.SpectrumIdentificationListType()
            {
                SpectrumIdentificationResult = new mzIdentML120.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML120.Generated.SpectrumIdentificationResultType()
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML120.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML120.Generated.SpectrumIdentificationItemType()
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
                {
                    new mzIdentML120.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML120.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML120.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML120.Generated.IonTypeType()
            {
                FragmentArray = new mzIdentML120.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML120.Generated.FragmentArrayType()
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML120.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML120.Generated.PeptideEvidenceRefType()
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.Inputs = new mzIdentML120.Generated.InputsType()
            {
                SpectraData = new mzIdentML120.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML120.Generated.SpectraDataType()
            {
                FileFormat = new mzIdentML120.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML120.Generated.CVParamType()
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML120.Generated.SequenceCollectionType()
            {
                PeptideEvidence = new mzIdentML120.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML120.Generated.PeptideEvidenceType()
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML120.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML120.Generated.PeptideType()
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML120.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML120.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML120.Generated.DBSequenceType()
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML120.Generated.ModificationType()
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML120.Generated.CVParamType()
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML120.Generated.AnalysisProtocolCollectionType()
            {
                SpectrumIdentificationProtocol = new mzIdentML120.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML120.Generated.SpectrumIdentificationProtocolType()
            {
                ParentTolerance = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML120.Generated.CVParamType()
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML120.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML120.Generated.CVParamType()
            {
                unitName = "dalton",
                value = "0.01"
            };

            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlWriter writer = XmlWriter.Create("myIdentifications.mzid", settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void MzidTest()
        {
            string version = "1.1.0";
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType110
            {
                version = version,
                id = "",
                DataCollection = new mzIdentML110.Generated.DataCollectionType()
            };
            _mzid.Provider = new mzIdentML110.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML110.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML110.Generated.RoleType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            _mzid.DataCollection.AnalysisData = new mzIdentML110.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML110.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML110.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML110.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML110.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML110.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML110.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
                {
                    new mzIdentML110.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML110.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML110.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML110.Generated.IonTypeType
            {
                FragmentArray = new mzIdentML110.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML110.Generated.FragmentArrayType
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML110.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML110.Generated.PeptideEvidenceRefType
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].passThreshold = true;

            _mzid.DataCollection.Inputs = new mzIdentML110.Generated.InputsType
            {
                SpectraData = new mzIdentML110.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML110.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML110.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML110.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML110.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML110.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML110.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                start = 2,
                end = 34,
                isDecoy = false,
                peptide_ref = "P_1",
                dBSequence_ref = "DB_1",
                id = "PE_1"
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML110.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML110.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML110.Generated.ModificationType[1],
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML110.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML110.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML110.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML110.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML110.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML110.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML110.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML110.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };

            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlWriter writer = XmlWriter.Create("myIdentifications.mzid", settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void Mzid110Test()
        {
            string version = "1.1.0";
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML110.Generated.MzIdentMLType110));
            var _mzid = new mzIdentML110.Generated.MzIdentMLType110
            {
                version = version,
                id = "",
                DataCollection = new mzIdentML110.Generated.DataCollectionType()
            };
            _mzid.Provider = new mzIdentML110.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML110.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML110.Generated.RoleType()
                    {
                        cvParam = new mzIdentML110.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };

            _mzid.DataCollection.AnalysisData = new mzIdentML110.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML110.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML110.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML110.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML110.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML110.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML110.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
                {
                    new mzIdentML110.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML110.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML110.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML110.Generated.IonTypeType
            {
                FragmentArray = new mzIdentML110.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML110.Generated.FragmentArrayType
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML110.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML110.Generated.PeptideEvidenceRefType
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.Inputs = new mzIdentML110.Generated.InputsType
            {
                SpectraData = new mzIdentML110.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML110.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML110.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML110.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML110.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML110.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML110.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML110.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML110.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML110.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML110.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML110.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML110.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML110.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML110.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML110.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML110.Generated.SpectrumIdentificationProtocolType()
            {
                ParentTolerance = new mzIdentML110.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML110.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML110.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };

            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlWriter writer = XmlWriter.Create("myIdentifications.mzid", settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void Mzid111Test_()
        {
            string version = "1.1.1";
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML111.Generated.MzIdentMLType111));
            var _mzid = new mzIdentML111.Generated.MzIdentMLType111
            {
                version = version,
                id = "",
                DataCollection = new mzIdentML111.Generated.DataCollectionType()
            };
            _mzid.Provider = new mzIdentML111.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML111.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML111.Generated.RoleType()
                    {
                        cvParam = new mzIdentML111.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML111.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML111.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML111.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML111.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML111.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML111.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML111.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML111.Generated.CVParamType[]
                {
                    new mzIdentML111.Generated.CVParamType
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] =
                new mzIdentML111.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation =
                new mzIdentML111.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] =
                new mzIdentML111.Generated.IonTypeType
                {
                    FragmentArray = new mzIdentML111.Generated.FragmentArrayType[1]
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] =
                new mzIdentML111.Generated.FragmentArrayType
                {
                    values = new float[3] { 200, 300, 400 }
                };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef =
                new mzIdentML111.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] =
                new mzIdentML111.Generated.PeptideEvidenceRefType
                {
                    peptideEvidence_ref = "PE_1"
                };
            _mzid.DataCollection.Inputs = new mzIdentML111.Generated.InputsType
            {
                SpectraData = new mzIdentML111.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML111.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML111.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML111.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML111.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML111.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML111.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML111.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML111.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML111.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML111.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML111.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML111.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML111.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML111.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML111.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML111.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML111.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML111.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML111.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };

            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlWriter writer = XmlWriter.Create("myIdentifications.mzid", settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }

        [Test]
        public void Mzid120Test_()
        {
            string version = "1.2.0";
            XmlSerializer _indexedSerializer = new XmlSerializer(typeof(mzIdentML120.Generated.MzIdentMLType120));
            var _mzid = new mzIdentML120.Generated.MzIdentMLType120
            {
                version = version,
                id = "",
                DataCollection = new mzIdentML120.Generated.DataCollectionType()
            };
            _mzid.Provider = new mzIdentML120.Generated.ProviderType()
            {
                id = "PROVIDER",
                ContactRole = new mzIdentML120.Generated.ContactRoleType()
                {
                    contact_ref = "UWMadisonSmithGroup",
                    Role = new mzIdentML120.Generated.RoleType()
                    {
                        cvParam = new mzIdentML120.Generated.CVParamType()
                        {
                            accession = "MS:1001271",
                            name = "researcher",
                            cvRef = "PSI-MS"
                        },
                    },
                },
            };
            _mzid.DataCollection.AnalysisData = new mzIdentML120.Generated.AnalysisDataType
            {
                SpectrumIdentificationList = new mzIdentML120.Generated.SpectrumIdentificationListType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0] = new mzIdentML120.Generated.SpectrumIdentificationListType
            {
                SpectrumIdentificationResult = new mzIdentML120.Generated.SpectrumIdentificationResultType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0] = new mzIdentML120.Generated.SpectrumIdentificationResultType
            {
                spectrumID = "spectrum 2",
                SpectrumIdentificationItem = new mzIdentML120.Generated.SpectrumIdentificationItemType[50]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0] = new mzIdentML120.Generated.SpectrumIdentificationItemType
            {
                experimentalMassToCharge = 1134.2609130203 + 0.000001 * 1134.2609130203 + 0.000001,
                calculatedMassToCharge = 1134.26091302033,
                calculatedMassToChargeSpecified = true,
                chargeState = 3,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
                {
                    new mzIdentML120.Generated.CVParamType()
                    {
                    accession = "MS:1002354",
                    value = "0.05"
                    }
                }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[1] = new mzIdentML120.Generated.SpectrumIdentificationItemType();
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation = new mzIdentML120.Generated.IonTypeType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0] = new mzIdentML120.Generated.IonTypeType
            {
                FragmentArray = new mzIdentML120.Generated.FragmentArrayType[1]
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].Fragmentation[0].FragmentArray[0] = new mzIdentML120.Generated.FragmentArrayType
            {
                values = new float[3] { 200, 300, 400 }
            };
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef = new mzIdentML120.Generated.PeptideEvidenceRefType[1];
            _mzid.DataCollection.AnalysisData.SpectrumIdentificationList[0].SpectrumIdentificationResult[0].SpectrumIdentificationItem[0].PeptideEvidenceRef[0] = new mzIdentML120.Generated.PeptideEvidenceRefType
            {
                peptideEvidence_ref = "PE_1"
            };
            _mzid.DataCollection.Inputs = new mzIdentML120.Generated.InputsType
            {
                SpectraData = new mzIdentML120.Generated.SpectraDataType[1]
            };
            _mzid.DataCollection.Inputs.SpectraData[0] = new mzIdentML120.Generated.SpectraDataType
            {
                FileFormat = new mzIdentML120.Generated.FileFormatType()
            };
            _mzid.DataCollection.Inputs.SpectraData[0].FileFormat.cvParam = new mzIdentML120.Generated.CVParamType
            {
                name = "mzML format"
            };
            _mzid.SequenceCollection = new mzIdentML120.Generated.SequenceCollectionType
            {
                PeptideEvidence = new mzIdentML120.Generated.PeptideEvidenceType[1]
            };
            _mzid.SequenceCollection.PeptideEvidence[0] = new mzIdentML120.Generated.PeptideEvidenceType
            {
                endSpecified = true,
                startSpecified = true,
                isDecoy = false,
                start = 2,
                end = 34,
                dBSequence_ref = "DB_1",
                peptide_ref = "P_1",
                id = "PE_1",
            };
            _mzid.SequenceCollection.Peptide = new mzIdentML120.Generated.PeptideType[1];
            _mzid.SequenceCollection.Peptide[0] = new mzIdentML120.Generated.PeptideType
            {
                id = "P_1",
                PeptideSequence = "GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR",
                Modification = new mzIdentML120.Generated.ModificationType[1]
            };
            _mzid.SequenceCollection.DBSequence = new mzIdentML120.Generated.DBSequenceType[1];
            _mzid.SequenceCollection.DBSequence[0] = new mzIdentML120.Generated.DBSequenceType
            {
                id = "DB_1",
                name = "Protein name",
                accession = "ACCESSION",
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0] = new mzIdentML120.Generated.ModificationType
            {
                locationSpecified = true,
                location = 17,
                monoisotopicMassDeltaSpecified = true,
                monoisotopicMassDelta = 57.02146373,
                cvParam = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.SequenceCollection.Peptide[0].Modification[0].cvParam[0] = new mzIdentML120.Generated.CVParamType
            {
                accession = "MS:1001460",
                name = "unknown modification",
                value = "Carbamidomethyl",
                cvRef = "PSI-MS"
            };
            _mzid.AnalysisProtocolCollection = new mzIdentML120.Generated.AnalysisProtocolCollectionType
            {
                SpectrumIdentificationProtocol = new mzIdentML120.Generated.SpectrumIdentificationProtocolType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0] = new mzIdentML120.Generated.SpectrumIdentificationProtocolType
            {
                ParentTolerance = new mzIdentML120.Generated.CVParamType[1]
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].ParentTolerance[0] = new mzIdentML120.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.1"
            };
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance = new mzIdentML120.Generated.CVParamType[1];
            _mzid.AnalysisProtocolCollection.SpectrumIdentificationProtocol[0].FragmentTolerance[0] = new mzIdentML120.Generated.CVParamType
            {
                unitName = "dalton",
                value = "0.01"
            };

            UTF8Encoding utf8EmitBOM = new UTF8Encoding(false);
            XmlWriterSettings settings = new XmlWriterSettings()
            {
                NewLineChars = "\n",
                Indent = true,
                Encoding = utf8EmitBOM,
            };
            XmlWriter writer = XmlWriter.Create("myIdentifications.mzid", settings);
            _indexedSerializer.Serialize(writer, _mzid);
            writer.Close();

            var identifications = new MzidIdentifications("myIdentifications.mzid");

            Assert.AreEqual(1134.26091302033, identifications.CalculatedMassToCharge(0, 0));
            Assert.AreEqual(3, identifications.ChargeState(0, 0));
            Assert.AreEqual(1, identifications.Count);
            Assert.AreEqual(1134.26091302033 + 0.000001 * 1134.2609130203 + 0.000001, identifications.ExperimentalMassToCharge(0, 0), 1e-10);
            Assert.IsFalse(identifications.IsDecoy(0, 0));
            Assert.AreEqual("MS:1001460", identifications.ModificationAcession(0, 0, 0));
            Assert.AreEqual("PSI-MS", identifications.ModificationDictionary(0, 0, 0));
            Assert.AreEqual("Carbamidomethyl", identifications.ModificationValue(0, 0, 0));
            Assert.AreEqual(17, identifications.ModificationLocation(0, 0, 0));
            Assert.AreEqual(57.02146373, identifications.ModificationMass(0, 0, 0));
            Assert.AreEqual("spectrum 2", identifications.Ms2SpectrumID(0));
            Assert.AreEqual(1, identifications.NumModifications(0, 0));
            Assert.AreEqual("GPEAPPPALPAGAPPPCTAVTSDHLNSLLGNILR", identifications.PeptideSequenceWithoutModifications(0, 0));
            Assert.AreEqual(0.1, identifications.ParentTolerance.Value);
            Assert.AreEqual(0.01, identifications.FragmentTolerance.Value);
            Assert.AreEqual(.05, identifications.QValue(0, 0));
            Assert.AreEqual("Protein name", identifications.ProteinFullName(0, 0));
            Assert.AreEqual("ACCESSION", identifications.ProteinAccession(0, 0));
            Assert.AreEqual(new float[3] { 200, 300, 400 }, identifications.MatchedIons(0, 0, 0));
            Assert.AreEqual(3, identifications.MatchedIonCounts(0, 0, 0));
            Assert.AreEqual("2", identifications.StartResidueInProtein(0, 0));
            Assert.AreEqual("34", identifications.EndResidueInProtein(0, 0));
            Assert.AreEqual(2, identifications.NumPSMsFromScan(0));
        }
    }
}