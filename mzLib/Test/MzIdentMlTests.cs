using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;
using mzIdentML111.Generated;
using System.Windows.Media;
using System.Xml.Serialization;
using System.IO;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class MzIdentMlTests
    {
        [Test]
        public void MzIdentMLType_Properties_SetAndGet()
        {
            var obj = new MzIdentMLType111();

            var cvList = new cvType[1] { new cvType() };
            var analysisSoftwareList = new AnalysisSoftwareType[1] { new AnalysisSoftwareType() };
            var provider = new ProviderType();
            var auditCollection = new AbstractContactType[1] { new OrganizationType() };
            var sampleCollection = new SampleType[1] { new SampleType() };
            var sequenceCollection = new SequenceCollectionType();
            var analysisCollection = new AnalysisCollectionType();
            var analysisProtocolCollection = new AnalysisProtocolCollectionType();
            var dataCollection = new DataCollectionType();
            var biblioRefs = new BibliographicReferenceType[1] { new BibliographicReferenceType() };
            var date = new DateTime(2022, 1, 1);

            obj.cvList = cvList;
            obj.AnalysisSoftwareList = analysisSoftwareList;
            obj.Provider = provider;
            obj.AuditCollection = auditCollection;
            obj.AnalysisSampleCollection = sampleCollection;
            obj.SequenceCollection = sequenceCollection;
            obj.AnalysisCollection = analysisCollection;
            obj.AnalysisProtocolCollection = analysisProtocolCollection;
            obj.DataCollection = dataCollection;
            obj.BibliographicReference = biblioRefs;
            obj.creationDate = date;
            obj.creationDateSpecified = true;
            obj.version = "1.2.3";

            Assert.That(obj.cvList, Is.EqualTo(cvList));
            Assert.That(obj.AnalysisSoftwareList, Is.EqualTo(analysisSoftwareList));
            Assert.That(obj.Provider, Is.EqualTo(provider));
            Assert.That(obj.AuditCollection, Is.EqualTo(auditCollection));
            Assert.That(obj.AnalysisSampleCollection, Is.EqualTo(sampleCollection));
            Assert.That(obj.SequenceCollection, Is.EqualTo(sequenceCollection));
            Assert.That(obj.AnalysisCollection, Is.EqualTo(analysisCollection));
            Assert.That(obj.AnalysisProtocolCollection, Is.EqualTo(analysisProtocolCollection));
            Assert.That(obj.DataCollection, Is.EqualTo(dataCollection));
            Assert.That(obj.BibliographicReference, Is.EqualTo(biblioRefs));
            Assert.That(obj.creationDate, Is.EqualTo(date));
            Assert.That(obj.creationDateSpecified, Is.EqualTo(true));
            Assert.That(obj.version, Is.EqualTo("1.2.3"));
        }

        [Test]
        public void cvType_Properties_SetAndGet()
        {
            var obj = new cvType();

            obj.fullName = "Full Name";
            obj.version = "2.0";
            obj.uri = "http://example.com";
            obj.id = "cv1";

            Assert.That(obj.fullName, Is.EqualTo("Full Name"));
            Assert.That(obj.version, Is.EqualTo("2.0"));
            Assert.That(obj.uri, Is.EqualTo("http://example.com"));
            Assert.That(obj.id, Is.EqualTo("cv1"));
        }

        [Test]
        public void AnalysisSoftwareType_Properties_SetAndGet()
        {
            var obj = new AnalysisSoftwareType();

            obj.name = "Software";
            obj.version = "1.0";
            obj.uri = "http://software.com";
            obj.id = "as1";

            Assert.That(obj.name, Is.EqualTo("Software"));
            Assert.That(obj.version, Is.EqualTo("1.0"));
            Assert.That(obj.uri, Is.EqualTo("http://software.com"));
            Assert.That(obj.id, Is.EqualTo("as1"));
        }

        [Test]
        public void ProviderType_Properties_SetAndGet()
        {
            var obj = new ProviderType();

            var contactRole = new ContactRoleType();
            obj.analysisSoftware_ref = "as1";
            obj.ContactRole = contactRole;

            Assert.That(obj.analysisSoftware_ref, Is.EqualTo("as1"));
            Assert.That(obj.ContactRole, Is.EqualTo(contactRole));
        }

        [Test]
        public void SampleType_Properties_SetAndGet()
        {
            var obj = new SampleType();

            obj.name = "Sample";
            obj.id = "sample1";
            obj.cvParam = new CVParamType[1] { new CVParamType() };
            obj.userParam = new UserParamType[1] { new UserParamType() };

            Assert.That(obj.name, Is.EqualTo("Sample"));
            Assert.That(obj.id, Is.EqualTo("sample1"));
            Assert.That(obj.cvParam, Is.EqualTo(new CVParamType[1] { obj.cvParam[0] }));
            Assert.That(obj.userParam, Is.EqualTo(new UserParamType[1] { obj.userParam[0] }));
        }

        [Test]
        public void SequenceCollectionType_Properties_SetAndGet()
        {
            var obj = new SequenceCollectionType();

            obj.DBSequence = new DBSequenceType[1] { new DBSequenceType() };
            obj.Peptide = new PeptideType[1] { new PeptideType() };
            obj.PeptideEvidence = new PeptideEvidenceType[1] { new PeptideEvidenceType() };

            Assert.That(obj.DBSequence, Is.EqualTo(new DBSequenceType[1] { obj.DBSequence[0] }));
            Assert.That(obj.Peptide, Is.EqualTo(new PeptideType[1] { obj.Peptide[0] }));
            Assert.That(obj.PeptideEvidence, Is.EqualTo(new PeptideEvidenceType[1] { obj.PeptideEvidence[0] }));
        }

        [Test]
        public void AnalysisCollectionType_Properties_SetAndGet()
        {
            var obj = new AnalysisCollectionType();

            obj.SpectrumIdentification = new SpectrumIdentificationType[1] { new SpectrumIdentificationType() };
            obj.ProteinDetection = new ProteinDetectionType();

            Assert.That(obj.SpectrumIdentification, Is.EqualTo(new SpectrumIdentificationType[1] { obj.SpectrumIdentification[0] }));
            Assert.That(obj.ProteinDetection, Is.EqualTo(obj.ProteinDetection));
        }

        [Test]
        public void AnalysisProtocolCollectionType_Properties_SetAndGet()
        {
            var obj = new AnalysisProtocolCollectionType();

            obj.SpectrumIdentificationProtocol = new SpectrumIdentificationProtocolType[1] { new SpectrumIdentificationProtocolType() };
            obj.ProteinDetectionProtocol = new ProteinDetectionProtocolType();

            Assert.That(obj.SpectrumIdentificationProtocol, Is.EqualTo(new SpectrumIdentificationProtocolType[1] { obj.SpectrumIdentificationProtocol[0] }));
            Assert.That(obj.ProteinDetectionProtocol, Is.EqualTo(obj.ProteinDetectionProtocol));
        }

        [Test]
        public void BibliographicReferenceType_Properties_SetAndGet()
        {
            var obj = new BibliographicReferenceType();

            obj.id = "b1";
            obj.authors = "Author";
            obj.title = "Title";
            obj.year = 2020;
            obj.volume = "1";
            obj.issue = "2";
            obj.pages = "3-4";
            obj.doi = "10.1000/xyz";
            obj.editor = "Editor";
            obj.publisher = "Publisher";

            Assert.That(obj.id, Is.EqualTo("b1"));
            Assert.That(obj.authors, Is.EqualTo("Author"));
            Assert.That(obj.title, Is.EqualTo("Title"));
            Assert.That(obj.year, Is.EqualTo(2020));
            Assert.That(obj.volume, Is.EqualTo("1"));
            Assert.That(obj.issue, Is.EqualTo("2"));
            Assert.That(obj.pages, Is.EqualTo("3-4"));
            Assert.That(obj.doi, Is.EqualTo("10.1000/xyz"));
            Assert.That(obj.editor, Is.EqualTo("Editor"));
            Assert.That(obj.publisher, Is.EqualTo("Publisher"));
        }
        [Test]
        public void ProteinDetectionHypothesisType_Properties_SetAndGet()
        {
            var obj = new ProteinDetectionHypothesisType();
            var peptideHypothesis = new PeptideHypothesisType[1] { new PeptideHypothesisType() };
            var cvParams = new CVParamType[1] { new CVParamType() };
            var userParams = new UserParamType[1] { new UserParamType() };

            obj.PeptideHypothesis = peptideHypothesis;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.dBSequence_ref = "dbseq";
            obj.passThreshold = true;

            Assert.That(obj.PeptideHypothesis, Is.EqualTo(peptideHypothesis));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.dBSequence_ref, Is.EqualTo("dbseq"));
            Assert.That(obj.passThreshold, Is.EqualTo(true));
        }

        [Test]
        public void ProteinAmbiguityGroupType_Properties_SetAndGet()
        {
            var obj = new ProteinAmbiguityGroupType();
            var pdh = new ProteinDetectionHypothesisType[1] { new ProteinDetectionHypothesisType() };
            var cvParams = new CVParamType[1] { new CVParamType() };
            var userParams = new UserParamType[1] { new UserParamType() };

            obj.ProteinDetectionHypothesis = pdh;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.ProteinDetectionHypothesis, Is.EqualTo(pdh));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void ProteinDetectionListType_Properties_SetAndGet()
        {
            var obj = new ProteinDetectionListType();
            var pag = new ProteinAmbiguityGroupType[1] { new ProteinAmbiguityGroupType() };
            var cvParams = new CVParamType[1] { new CVParamType() };
            var userParams = new UserParamType[1] { new UserParamType() };

            obj.ProteinAmbiguityGroup = pag;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.ProteinAmbiguityGroup, Is.EqualTo(pag));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void SpectrumIdentificationItemType_Properties_SetAndGet()
        {
            var obj = new SpectrumIdentificationItemType();
            var peptideEvidenceRefs = new PeptideEvidenceRefType[1] { new PeptideEvidenceRefType() };
            var fragmentation = new IonTypeType[1] { new IonTypeType() };
            var cvParams = new CVParamType[1] { new CVParamType() };
            var userParams = new UserParamType[1] { new UserParamType() };

            obj.PeptideEvidenceRef = peptideEvidenceRefs;
            obj.Fragmentation = fragmentation;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.chargeState = 2;
            obj.experimentalMassToCharge = 123.45;
            obj.calculatedMassToCharge = 120.12;
            obj.calculatedMassToChargeSpecified = true;
            obj.calculatedPI = 6.5f;
            obj.calculatedPISpecified = true;
            obj.peptide_ref = "pep1";
            obj.rank = 1;
            obj.passThreshold = true;
            obj.massTable_ref = "mt1";
            obj.sample_ref = "s1";

            Assert.That(obj.PeptideEvidenceRef, Is.EqualTo(peptideEvidenceRefs));
            Assert.That(obj.Fragmentation, Is.EqualTo(fragmentation));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.chargeState, Is.EqualTo(2));
            Assert.That(obj.experimentalMassToCharge, Is.EqualTo(123.45));
            Assert.That(obj.calculatedMassToCharge, Is.EqualTo(120.12));
            Assert.That(obj.calculatedMassToChargeSpecified, Is.EqualTo(true));
            Assert.That(obj.calculatedPI, Is.EqualTo(6.5f));
            Assert.That(obj.calculatedPISpecified, Is.EqualTo(true));
            Assert.That(obj.peptide_ref, Is.EqualTo("pep1"));
            Assert.That(obj.rank, Is.EqualTo(1));
            Assert.That(obj.passThreshold, Is.EqualTo(true));
            Assert.That(obj.massTable_ref, Is.EqualTo("mt1"));
            Assert.That(obj.sample_ref, Is.EqualTo("s1"));
        }

        [Test]
        public void SpectrumIdentificationResultType_Properties_SetAndGet()
        {
            var obj = new SpectrumIdentificationResultType();
            var sii = new SpectrumIdentificationItemType[1] { new SpectrumIdentificationItemType() };
            var cvParams = new CVParamType[1] { new CVParamType() };
            var userParams = new UserParamType[1] { new UserParamType() };

            obj.SpectrumIdentificationItem = sii;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.spectrumID = "spec1";
            obj.spectraData_ref = "sd1";

            Assert.That(obj.SpectrumIdentificationItem, Is.EqualTo(sii));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.spectrumID, Is.EqualTo("spec1"));
            Assert.That(obj.spectraData_ref, Is.EqualTo("sd1"));
        }

        [Test]
        public void ExternalDataType_Properties_SetAndGet()
        {
            var obj = new ExternalDataType();
            var fileFormat = new FileFormatType();

            obj.ExternalFormatDocumentation = "http://doc";
            obj.FileFormat = fileFormat;
            obj.location = "http://location";

            Assert.That(obj.ExternalFormatDocumentation, Is.EqualTo("http://doc"));
            Assert.That(obj.FileFormat, Is.EqualTo(fileFormat));
            Assert.That(obj.location, Is.EqualTo("http://location"));
        }

        [Test]
        public void FileFormatType_Properties_SetAndGet()
        {
            var obj = new FileFormatType();
            var cvParam = new CVParamType();

            obj.cvParam = cvParam;

            Assert.That(obj.cvParam, Is.EqualTo(cvParam));
        }

        [Test]
        public void SpectraDataType_Properties_SetAndGet()
        {
            var obj = new SpectraDataType();
            var spectrumIDFormat = new SpectrumIDFormatType();

            obj.SpectrumIDFormat = spectrumIDFormat;

            Assert.That(obj.SpectrumIDFormat, Is.EqualTo(spectrumIDFormat));
        }

        [Test]
        public void SpectrumIDFormatType_Properties_SetAndGet()
        {
            var obj = new SpectrumIDFormatType();
            var cvParam = new CVParamType();

            obj.cvParam = cvParam;

            Assert.That(obj.cvParam, Is.EqualTo(cvParam));
        }

        [Test]
        public void SourceFileType_Properties_SetAndGet()
        {
            var obj = new SourceFileType();
            var cvParams = new CVParamType[1] { new CVParamType() };
            var userParams = new UserParamType[1] { new UserParamType() };

            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }
        [Test]
        public void PeptideEvidenceRefType_Properties_SetAndGet()
        {
            var obj = new PeptideEvidenceRefType();
            obj.peptideEvidence_ref = "pe1";
            Assert.That(obj.peptideEvidence_ref, Is.EqualTo("pe1"));
        }

        [Test]
        public void CVParamType_Properties_SetAndGet()
        {
            var obj = new CVParamType();
            obj.cvRef = "cv1";
            obj.accession = "ACC";
            obj.value = "val";
            obj.unitCvRef = "unitcv";
            obj.unitAccession = "unitacc";
            obj.unitName = "unitname";
            obj.name = "paramname";

            Assert.That(obj.cvRef, Is.EqualTo("cv1"));
            Assert.That(obj.accession, Is.EqualTo("ACC"));
            Assert.That(obj.value, Is.EqualTo("val"));
            Assert.That(obj.unitCvRef, Is.EqualTo("unitcv"));
            Assert.That(obj.unitAccession, Is.EqualTo("unitacc"));
            Assert.That(obj.unitName, Is.EqualTo("unitname"));
            Assert.That(obj.name, Is.EqualTo("paramname"));
        }

        [Test]
        public void UserParamType_Properties_SetAndGet()
        {
            var obj = new UserParamType();
            obj.name = "username";
            obj.value = "uservalue";
            obj.type = "usertype";
            obj.unitAccession = "unitacc";
            obj.unitCvRef = "unitcv";
            obj.unitName = "unitname";

            Assert.That(obj.name, Is.EqualTo("username"));
            Assert.That(obj.value, Is.EqualTo("uservalue"));
            Assert.That(obj.type, Is.EqualTo("usertype"));
            Assert.That(obj.unitAccession, Is.EqualTo("unitacc"));
            Assert.That(obj.unitCvRef, Is.EqualTo("unitcv"));
            Assert.That(obj.unitName, Is.EqualTo("unitname"));
        }

        [Test]
        public void PeptideEvidenceType_Properties_SetAndGet()
        {
            var obj = new PeptideEvidenceType();
            obj.id = "pe1";
            obj.peptide_ref = "pep1";
            obj.dBSequence_ref = "db1";
            obj.start = 1;
            obj.end = 10;
            obj.pre = "K";
            obj.post = "R";
            obj.frame = 2;
            obj.isDecoy = true;
            obj.translationTable_ref = "tt1";
            obj.name = "evidence";

            Assert.That(obj.id, Is.EqualTo("pe1"));
            Assert.That(obj.peptide_ref, Is.EqualTo("pep1"));
            Assert.That(obj.dBSequence_ref, Is.EqualTo("db1"));
            Assert.That(obj.start, Is.EqualTo(1));
            Assert.That(obj.end, Is.EqualTo(10));
            Assert.That(obj.pre, Is.EqualTo("K"));
            Assert.That(obj.post, Is.EqualTo("R"));
            Assert.That(obj.frame, Is.EqualTo(2));
            Assert.That(obj.isDecoy, Is.EqualTo(true));
            Assert.That(obj.translationTable_ref, Is.EqualTo("tt1"));
            Assert.That(obj.name, Is.EqualTo("evidence"));
        }

        [Test]
        public void DBSequenceType_Properties_SetAndGet()
        {
            var obj = new DBSequenceType();
            obj.id = "db1";
            obj.length = 100;
            obj.lengthSpecified = true;
            obj.searchDatabase_ref = "sdb1";
            obj.accession = "acc1";
            obj.name = "dbseq";
            obj.cvParam = new CVParamType[1] { new CVParamType() };
            obj.userParam = new UserParamType[1] { new UserParamType() };

            Assert.That(obj.id, Is.EqualTo("db1"));
            Assert.That(obj.length, Is.EqualTo(100));
            Assert.That(obj.lengthSpecified, Is.EqualTo(true));
            Assert.That(obj.searchDatabase_ref, Is.EqualTo("sdb1"));
            Assert.That(obj.accession, Is.EqualTo("acc1"));
            Assert.That(obj.name, Is.EqualTo("dbseq"));
            Assert.That(obj.cvParam, Is.EqualTo(new CVParamType[1] { obj.cvParam[0] }));
            Assert.That(obj.userParam, Is.EqualTo(new UserParamType[1] { obj.userParam[0] }));
        }

        [Test]
        public void PeptideType_Properties_SetAndGet()
        {
            var obj = new PeptideType();
            obj.id = "pep1";
            obj.name = "peptide";
            obj.cvParam = new CVParamType[1] { new CVParamType() };
            obj.userParam = new UserParamType[1] { new UserParamType() };

            Assert.That(obj.id, Is.EqualTo("pep1"));
            Assert.That(obj.name, Is.EqualTo("peptide"));
            Assert.That(obj.cvParam, Is.EqualTo(new CVParamType[1] { obj.cvParam[0] }));
            Assert.That(obj.userParam, Is.EqualTo(new UserParamType[1] { obj.userParam[0] }));
        }

        [Test]
        public void ModificationType_Properties_SetAndGet()
        {
            var obj = new ModificationType();
            obj.location = 3;
            obj.locationSpecified = true;
            obj.residues = new string[2] { "A", "C" };
            obj.cvParam = new CVParamType[1] { new CVParamType() };

            Assert.That(obj.location, Is.EqualTo(3));
            Assert.That(obj.locationSpecified, Is.EqualTo(true));
            Assert.That(obj.residues, Is.EqualTo(new string[2] { "A", "C" }));
            Assert.That(obj.cvParam, Is.EqualTo(new CVParamType[1] { obj.cvParam[0] }));
        }
        [Test]
        public void SpectrumIdentificationItemType_Property_Set_Get_Works()
        {
            var item = new SpectrumIdentificationItemType();

            var peptideEvidenceRefs = new PeptideEvidenceRefType[] { new PeptideEvidenceRefType { peptideEvidence_ref = "pep1" } };
            var fragmentations = new IonTypeType[] { new IonTypeType { charge = 2 } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            item.PeptideEvidenceRef = peptideEvidenceRefs;
            item.Fragmentation = fragmentations;
            item.cvParam = cvParams;
            item.userParam = userParams;
            item.chargeState = 3;
            item.experimentalMassToCharge = 123.45;
            item.calculatedMassToCharge = 120.12;
            item.calculatedMassToChargeSpecified = true;
            item.calculatedPI = 6.5f;
            item.calculatedPISpecified = true;
            item.peptide_ref = "peptide1";
            item.rank = 1;
            item.passThreshold = true;
            item.massTable_ref = "mt1";
            item.sample_ref = "sample1";

            Assert.That(item.PeptideEvidenceRef, Is.EqualTo(peptideEvidenceRefs));
            Assert.That(item.Fragmentation, Is.EqualTo(fragmentations));
            Assert.That(item.cvParam, Is.EqualTo(cvParams));
            Assert.That(item.userParam, Is.EqualTo(userParams));
            Assert.That(item.chargeState, Is.EqualTo(3));
            Assert.That(item.experimentalMassToCharge, Is.EqualTo(123.45));
            Assert.That(item.calculatedMassToCharge, Is.EqualTo(120.12));
            Assert.That(item.calculatedMassToChargeSpecified, Is.EqualTo(true));
            Assert.That(item.calculatedPI, Is.EqualTo(6.5f));
            Assert.That(item.calculatedPISpecified, Is.EqualTo(true));
            Assert.That(item.peptide_ref, Is.EqualTo("peptide1"));
            Assert.That(item.rank, Is.EqualTo(1));
            Assert.That(item.passThreshold, Is.EqualTo(true));
            Assert.That(item.massTable_ref, Is.EqualTo("mt1"));
            Assert.That(item.sample_ref, Is.EqualTo("sample1"));
        }

        [Test]
        public void SpectrumIdentificationResultType_Property_Set_Get_Works()
        {
            var result = new SpectrumIdentificationResultType();

            var items = new SpectrumIdentificationItemType[] { new SpectrumIdentificationItemType { rank = 1 } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC2" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type2" } };

            result.SpectrumIdentificationItem = items;
            result.cvParam = cvParams;
            result.userParam = userParams;
            result.spectrumID = "spec1";
            result.spectraData_ref = "data1";

            Assert.That(result.SpectrumIdentificationItem, Is.EqualTo(items));
            Assert.That(result.cvParam, Is.EqualTo(cvParams));
            Assert.That(result.userParam, Is.EqualTo(userParams));
            Assert.That(result.spectrumID, Is.EqualTo("spec1"));
            Assert.That(result.spectraData_ref, Is.EqualTo("data1"));
        }
        [Test]
        public void ProteinDetectionHypothesisType_Properties_SetAndGetAlternate()
        {
            var obj = new ProteinDetectionHypothesisType();
            var peptideHypotheses = new PeptideHypothesisType[] { new PeptideHypothesisType { peptideEvidence_ref = "pe1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            obj.PeptideHypothesis = peptideHypotheses;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.id = "pdh1";
            obj.dBSequence_ref = "dbseq1";
            obj.passThreshold = true;

            Assert.That(obj.PeptideHypothesis, Is.EqualTo(peptideHypotheses));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.id, Is.EqualTo("pdh1"));
            Assert.That(obj.dBSequence_ref, Is.EqualTo("dbseq1"));
            Assert.That(obj.passThreshold, Is.EqualTo(true));
        }

        [Test]
        public void PeptideHypothesisType_Properties_SetAndGet()
        {
            var obj = new PeptideHypothesisType();
            var spectrumRefs = new SpectrumIdentificationItemRefType[] { new SpectrumIdentificationItemRefType { spectrumIdentificationItem_ref = "sii1" } };

            obj.SpectrumIdentificationItemRef = spectrumRefs;
            obj.peptideEvidence_ref = "pe1";

            Assert.That(obj.SpectrumIdentificationItemRef, Is.EqualTo(spectrumRefs));
            Assert.That(obj.peptideEvidence_ref, Is.EqualTo("pe1"));
        }

        [Test]
        public void SpectrumIdentificationItemRefType_Properties_SetAndGet()
        {
            var obj = new SpectrumIdentificationItemRefType();
            obj.spectrumIdentificationItem_ref = "sii1";

            Assert.That(obj.spectrumIdentificationItem_ref, Is.EqualTo("sii1"));
        }

        [Test]
        public void ProteinAmbiguityGroupType_Properties_SetAndGetAlternate()
        {
            var obj = new ProteinAmbiguityGroupType();
            var pdh = new ProteinDetectionHypothesisType[] { new ProteinDetectionHypothesisType { id = "pdh1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC2" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type2" } };

            obj.ProteinDetectionHypothesis = pdh;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.id = "pag1";

            Assert.That(obj.ProteinDetectionHypothesis, Is.EqualTo(pdh));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.id, Is.EqualTo("pag1"));
        }

        [Test]
        public void ProteinDetectionListType_Properties_SetAndGetAlternate()
        {
            var obj = new ProteinDetectionListType();
            var pag = new ProteinAmbiguityGroupType[] { new ProteinAmbiguityGroupType { id = "pag1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC3" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type3" } };

            obj.ProteinAmbiguityGroup = pag;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.id = "pdl1";

            Assert.That(obj.ProteinAmbiguityGroup, Is.EqualTo(pag));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.id, Is.EqualTo("pdl1"));
        }

        [Test]
        public void MzIdentMLType_Properties_SetAndGetAlternate()
        {
            var obj = new MzIdentMLType111();
            obj.id = "mzid1";
            obj.version = "1.1.1";
            obj.creationDate = new System.DateTime(2024, 6, 1);
            obj.creationDateSpecified = true;

            Assert.That(obj.id, Is.EqualTo("mzid1"));
            Assert.That(obj.version, Is.EqualTo("1.1.1"));
            Assert.That(obj.creationDate, Is.EqualTo(new System.DateTime(2024, 6, 1)));
            Assert.That(obj.creationDateSpecified, Is.EqualTo(true));
        }
        [Test]
        public void AnalysisProtocolCollectionType_Properties_SetAndGetA()
        {
            var obj = new AnalysisProtocolCollectionType();
            var spectrumIdentProtocols = new SpectrumIdentificationProtocolType[] { new SpectrumIdentificationProtocolType { id = "sip1" } };
            var proteinDetectionProtocol = new ProteinDetectionProtocolType { id = "pdp1" };

            obj.SpectrumIdentificationProtocol = spectrumIdentProtocols;
            obj.ProteinDetectionProtocol = proteinDetectionProtocol;

            Assert.That(obj.SpectrumIdentificationProtocol, Is.EqualTo(spectrumIdentProtocols));
            Assert.That(obj.ProteinDetectionProtocol, Is.EqualTo(proteinDetectionProtocol));
        }


        [Test]
        public void ProteinDetectionProtocolType_Properties_SetAndGet()
        {
            var obj = new ProteinDetectionProtocolType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC2" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type2" } };

            obj.id = "pdp1";
            obj.analysisSoftware_ref = "asw1";

            Assert.That(obj.id, Is.EqualTo("pdp1"));
            Assert.That(obj.analysisSoftware_ref, Is.EqualTo("asw1"));
        }

        [Test]
        public void SearchTypeType_Properties_SetAndGet()
        {
            var obj = new SearchDatabaseType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC3" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type3" } };

            obj.cvParam = cvParams;

            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
        }

        [Test]
        public void EnzymesType_Properties_SetAndGet()
        {
            var obj = new EnzymesType();
            var enzymes = new EnzymeType[] { new EnzymeType { id = "enz2" } };

            obj.Enzyme = enzymes;

            Assert.That(obj.Enzyme, Is.EqualTo(enzymes));
        }

        [Test]
        public void EnzymeType_Properties_SetAndGet()
        {
            var obj = new EnzymeType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC5" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type5" } };

            obj.id = "enz3";
            obj.name = "Trypsin";
            obj.nTermGain = "H";
            obj.cTermGain = "OH";
            obj.minDistance = 1;
            obj.minDistanceSpecified = true;
            obj.semiSpecific = true;
            obj.semiSpecificSpecified = true;

            Assert.That(obj.id, Is.EqualTo("enz3"));
            Assert.That(obj.name, Is.EqualTo("Trypsin"));
            Assert.That(obj.nTermGain, Is.EqualTo("H"));
            Assert.That(obj.cTermGain, Is.EqualTo("OH"));
            Assert.That(obj.minDistance, Is.EqualTo(1));
            Assert.That(obj.minDistanceSpecified, Is.EqualTo(true));
            Assert.That(obj.semiSpecific, Is.EqualTo(true));
            Assert.That(obj.semiSpecificSpecified, Is.EqualTo(true));
        }
        [Test]
        public void ModificationType_Properties_SetAndGetA()
        {
            var obj = new ModificationType();
            var residues = new string[] { "M" };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            obj.location = 3;
            obj.locationSpecified = true;
            obj.residues = residues;
            obj.avgMassDelta = 15.9949;
            obj.avgMassDeltaSpecified = true;
            obj.monoisotopicMassDelta = 15.9949;
            obj.monoisotopicMassDeltaSpecified = true;
            obj.cvParam = cvParams;

            Assert.That(obj.location, Is.EqualTo(3));
            Assert.That(obj.locationSpecified, Is.EqualTo(true));
            Assert.That(obj.residues, Is.EqualTo(residues));
            Assert.That(obj.avgMassDelta, Is.EqualTo(15.9949));
            Assert.That(obj.avgMassDeltaSpecified, Is.EqualTo(true));
            Assert.That(obj.monoisotopicMassDelta, Is.EqualTo(15.9949));
            Assert.That(obj.monoisotopicMassDeltaSpecified, Is.EqualTo(true));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
        }

        [Test]
        public void CVParamType_Properties_SetAndGetA()
        {
            var obj = new CVParamType();

            obj.cvRef = "MS";
            obj.accession = "ACC2";
            obj.value = "value";
            obj.unitCvRef = "UO";
            obj.unitAccession = "UO:0000000";
            obj.unitName = "unit";

            Assert.That(obj.cvRef, Is.EqualTo("MS"));
            Assert.That(obj.accession, Is.EqualTo("ACC2"));
            Assert.That(obj.value, Is.EqualTo("value"));
            Assert.That(obj.unitCvRef, Is.EqualTo("UO"));
            Assert.That(obj.unitAccession, Is.EqualTo("UO:0000000"));
            Assert.That(obj.unitName, Is.EqualTo("unit"));
        }

        [Test]
        public void UserParamType_Properties_SetAndGetA()
        {
            var obj = new UserParamType();

            obj.name = "param";
            obj.type = "string";
            obj.value = "test";
            obj.unitAccession = "UO:0000001";
            obj.unitCvRef = "UO";
            obj.unitName = "unit";

            Assert.That(obj.name, Is.EqualTo("param"));
            Assert.That(obj.type, Is.EqualTo("string"));
            Assert.That(obj.value, Is.EqualTo("test"));
            Assert.That(obj.unitAccession, Is.EqualTo("UO:0000001"));
            Assert.That(obj.unitCvRef, Is.EqualTo("UO"));
            Assert.That(obj.unitName, Is.EqualTo("unit"));
        }
        [Test]
        public void PeptideType_Properties_SetAndGetA()
        {
            var obj = new PeptideType();
            var modifications = new ModificationType[] { new ModificationType { location = 1 } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            obj.id = "pep1";
            obj.PeptideSequence = "ACDEFGHIK";
            obj.Modification = modifications;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("pep1"));
            Assert.That(obj.PeptideSequence, Is.EqualTo("ACDEFGHIK"));
            Assert.That(obj.Modification, Is.EqualTo(modifications));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void PeptideEvidenceType_Properties_SetAndGetA()
        {
            var obj = new PeptideEvidenceType();

            obj.id = "pe1";
            obj.peptide_ref = "pep1";
            obj.dBSequence_ref = "db1";
            obj.start = 2;
            obj.end = 10;
            obj.pre = "K";
            obj.post = "R";
            obj.frame = 1;
            obj.frameSpecified = true;
            obj.translationTable_ref = "tt1";
            obj.isDecoy = true;

            Assert.That(obj.id, Is.EqualTo("pe1"));
            Assert.That(obj.peptide_ref, Is.EqualTo("pep1"));
            Assert.That(obj.dBSequence_ref, Is.EqualTo("db1"));
            Assert.That(obj.start, Is.EqualTo(2));
            Assert.That(obj.end, Is.EqualTo(10));
            Assert.That(obj.pre, Is.EqualTo("K"));
            Assert.That(obj.post, Is.EqualTo("R"));
            Assert.That(obj.frame, Is.EqualTo(1));
            Assert.That(obj.frameSpecified, Is.EqualTo(true));
            Assert.That(obj.translationTable_ref, Is.EqualTo("tt1"));
            Assert.That(obj.isDecoy, Is.EqualTo(true));
        }

        [Test]
        public void DBSequenceType_Properties_SetAndGetA()
        {
            var obj = new DBSequenceType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC2" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type2" } };

            obj.id = "db1";
            obj.accession = "P12345";
            obj.searchDatabase_ref = "sdb1";
            obj.length = 100;
            obj.lengthSpecified = true;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("db1"));
            Assert.That(obj.accession, Is.EqualTo("P12345"));
            Assert.That(obj.searchDatabase_ref, Is.EqualTo("sdb1"));
            Assert.That(obj.length, Is.EqualTo(100));
            Assert.That(obj.lengthSpecified, Is.EqualTo(true));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void SpectraDataType_Properties_SetAndGetA()
        {
            var obj = new SpectraDataType();
            var userParams = new UserParamType[] { new UserParamType { type = "type3" } };
            var fileFormat = new FileFormatType { cvParam = new CVParamType { accession = "ACC3" } };
            var spectrumIDFormat = new SpectrumIDFormatType { cvParam = new CVParamType { accession = "ACC3" } };

            obj.id = "sd1";
            obj.location = "file://spectra.mzML";
            obj.name = "TestSpectra";
            obj.FileFormat = fileFormat;
            obj.SpectrumIDFormat = spectrumIDFormat;

            Assert.That(obj.id, Is.EqualTo("sd1"));
            Assert.That(obj.location, Is.EqualTo("file://spectra.mzML"));
            Assert.That(obj.name, Is.EqualTo("TestSpectra"));
            Assert.That(obj.FileFormat, Is.EqualTo(fileFormat));
            Assert.That(obj.SpectrumIDFormat, Is.EqualTo(spectrumIDFormat));
        }

        [Test]
        public void SourceFileType_Properties_SetAndGetA()
        {
            var obj = new SourceFileType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            obj.id = "sf1";
            obj.location = "file://source.raw";
            obj.name = "SourceFile";
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("sf1"));
            Assert.That(obj.location, Is.EqualTo("file://source.raw"));
            Assert.That(obj.name, Is.EqualTo("SourceFile"));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void PersonType_Properties_SetAndGet()
        {
            var obj = new PersonType();
            var affiliations = new AffiliationType[] { new AffiliationType { organization_ref = "org1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC4" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type4" } };

            obj.id = "person1";
            obj.name = "John Doe";
            obj.Affiliation = affiliations;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("person1"));
            Assert.That(obj.name, Is.EqualTo("John Doe"));
            Assert.That(obj.Affiliation, Is.EqualTo(affiliations));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void AffiliationType_Properties_SetAndGet()
        {
            var obj = new AffiliationType();

            obj.organization_ref = "org1";

            Assert.That(obj.organization_ref, Is.EqualTo("org1"));
        }
        [Test]
        public void OrganizationType_Properties_SetAndGet()
        {
            var obj = new OrganizationType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            obj.id = "org1";
            obj.name = "Institute of Science";
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("org1"));
            Assert.That(obj.name, Is.EqualTo("Institute of Science"));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void SampleType_Properties_SetAndGetA()
        {
            var obj = new SampleType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC2" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type2" } };

            obj.id = "sample1";
            obj.name = "Sample A";
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("sample1"));
            Assert.That(obj.name, Is.EqualTo("Sample A"));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }
        [Test]
        public void AnalysisDataType_Properties_SetAndGet()
        {
            var obj = new AnalysisDataType();
            var spectrumIdentList = new SpectrumIdentificationListType[] { new SpectrumIdentificationListType { id = "sil1" } };
            var proteinDetectionList = new ProteinDetectionListType { id = "pdl1" };

            obj.SpectrumIdentificationList = spectrumIdentList;
            obj.ProteinDetectionList = proteinDetectionList;

            Assert.That(obj.SpectrumIdentificationList, Is.EqualTo(spectrumIdentList));
            Assert.That(obj.ProteinDetectionList, Is.EqualTo(proteinDetectionList));
        }

        [Test]
        public void SpectrumIdentificationListType_Properties_SetAndGet()
        {
            var obj = new SpectrumIdentificationListType();
            var spectrumIdentResults = new SpectrumIdentificationResultType[] { new SpectrumIdentificationResultType { id = "sir1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            obj.id = "sil1";
            obj.numSequencesSearched = 1000;
            obj.numSequencesSearchedSpecified = true;
            obj.SpectrumIdentificationResult = spectrumIdentResults;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("sil1"));
            Assert.That(obj.numSequencesSearched, Is.EqualTo(1000));
            Assert.That(obj.numSequencesSearchedSpecified, Is.EqualTo(true));
            Assert.That(obj.SpectrumIdentificationResult, Is.EqualTo(spectrumIdentResults));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void SpectrumIdentificationResultType_Properties_SetAndGetA()
        {
            var obj = new SpectrumIdentificationResultType();
            var spectrumIdentItems = new SpectrumIdentificationItemType[] { new SpectrumIdentificationItemType { id = "sii1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC2" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type2" } };

            obj.id = "sir1";
            obj.spectrumID = "spec1";
            obj.spectraData_ref = "sd1";
            obj.SpectrumIdentificationItem = spectrumIdentItems;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("sir1"));
            Assert.That(obj.spectrumID, Is.EqualTo("spec1"));
            Assert.That(obj.spectraData_ref, Is.EqualTo("sd1"));
            Assert.That(obj.SpectrumIdentificationItem, Is.EqualTo(spectrumIdentItems));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void SpectrumIdentificationItemType_Properties_SetAndGetA()
        {
            var obj = new SpectrumIdentificationItemType();
            var peptideEvidenceRefs = new PeptideEvidenceRefType[] { new PeptideEvidenceRefType { peptideEvidence_ref = "pe1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC3" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type3" } };

            obj.id = "sii1";
            obj.chargeState = 2;
            obj.experimentalMassToCharge = 500.2;
            obj.calculatedMassToCharge = 500.1;
            obj.peptide_ref = "pep1";
            obj.rank = 1;
            obj.passThreshold = true;
            obj.PeptideEvidenceRef = peptideEvidenceRefs;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("sii1"));
            Assert.That(obj.chargeState, Is.EqualTo(2));
            Assert.That(obj.experimentalMassToCharge, Is.EqualTo(500.2));
            Assert.That(obj.calculatedMassToCharge, Is.EqualTo(500.1));
            Assert.That(obj.peptide_ref, Is.EqualTo("pep1"));
            Assert.That(obj.rank, Is.EqualTo(1));
            Assert.That(obj.passThreshold, Is.EqualTo(true));
            Assert.That(obj.PeptideEvidenceRef, Is.EqualTo(peptideEvidenceRefs));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void PeptideEvidenceRefType_Properties_SetAndGetA()
        {
            var obj = new PeptideEvidenceRefType();

            obj.peptideEvidence_ref = "pe1";

            Assert.That(obj.peptideEvidence_ref, Is.EqualTo("pe1"));
        }

        [Test]
        public void ProteinDetectionListType_Properties_SetAndGetA()
        {
            var obj = new ProteinDetectionListType();
            var proteinAmbiguityGroups = new ProteinAmbiguityGroupType[] { new ProteinAmbiguityGroupType { id = "pag1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC4" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type4" } };

            obj.id = "pdl1";
            obj.ProteinAmbiguityGroup = proteinAmbiguityGroups;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("pdl1"));
            Assert.That(obj.ProteinAmbiguityGroup, Is.EqualTo(proteinAmbiguityGroups));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void ProteinAmbiguityGroupType_Properties_SetAndGetA()
        {
            var obj = new ProteinAmbiguityGroupType();
            var proteinDetectionHypotheses = new ProteinDetectionHypothesisType[] { new ProteinDetectionHypothesisType { id = "pdh1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC5" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type5" } };

            obj.id = "pag1";
            obj.ProteinDetectionHypothesis = proteinDetectionHypotheses;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("pag1"));
            Assert.That(obj.ProteinDetectionHypothesis, Is.EqualTo(proteinDetectionHypotheses));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void ProteinDetectionHypothesisType_Properties_SetAndGetA()
        {
            var obj = new ProteinDetectionHypothesisType();
            var peptideHypothesis = new PeptideHypothesisType[] { new PeptideHypothesisType { peptideEvidence_ref = "pe1" } };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC6" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type6" } };

            obj.id = "pdh1";
            obj.dBSequence_ref = "db1";
            obj.PeptideHypothesis = peptideHypothesis;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("pdh1"));
            Assert.That(obj.dBSequence_ref, Is.EqualTo("db1"));
            Assert.That(obj.PeptideHypothesis, Is.EqualTo(peptideHypothesis));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void AnalysisProtocolCollectionType_Properties_SetAndGetB()
        {
            var obj = new AnalysisProtocolCollectionType();
            var spectrumIdentificationProtocol = new SpectrumIdentificationProtocolType[] { new SpectrumIdentificationProtocolType { id = "sip1" } };
            var proteinDetectionProtocol = new ProteinDetectionProtocolType { id = "pdp1" };

            obj.SpectrumIdentificationProtocol = spectrumIdentificationProtocol;
            obj.ProteinDetectionProtocol = proteinDetectionProtocol;

            Assert.That(obj.SpectrumIdentificationProtocol, Is.EqualTo(spectrumIdentificationProtocol));
            Assert.That(obj.ProteinDetectionProtocol, Is.EqualTo(proteinDetectionProtocol));
        }


        [Test]
        public void EnzymesType_Properties_SetAndGetA()
        {
            var obj = new EnzymesType();
            var enzymes = new EnzymeType[] { new EnzymeType { id = "enz2" } };

            obj.Enzyme = enzymes;

            Assert.That(obj.Enzyme, Is.EqualTo(enzymes));
        }

        [Test]
        public void EnzymeType_Properties_SetAndGetA()
        {
            var obj = new EnzymeType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC4" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type4" } };

            obj.id = "enz1";
            obj.name = "Trypsin";
            obj.SiteRegexp = "[KR]|[X]";
            obj.nTermGain = "H";
            obj.cTermGain = "OH";
            obj.minDistance = 1;
            obj.minDistanceSpecified = true;

            Assert.That(obj.id, Is.EqualTo("enz1"));
            Assert.That(obj.name, Is.EqualTo("Trypsin"));
            Assert.That(obj.SiteRegexp, Is.EqualTo("[KR]|[X]"));
            Assert.That(obj.nTermGain, Is.EqualTo("H"));
            Assert.That(obj.cTermGain, Is.EqualTo("OH"));
            Assert.That(obj.minDistance, Is.EqualTo(1));
            Assert.That(obj.minDistanceSpecified, Is.EqualTo(true));
        }

        [Test]
        public void MassTableType_Properties_SetAndGet()
        {
            var obj = new MassTableType();
            var residues = new string[] { "A", "C" };
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC6" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type6" } };

            obj.id = "mt1";
            obj.name = "MassTable1";
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("mt1"));
            Assert.That(obj.name, Is.EqualTo("MassTable1"));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }


        [Test]
        public void ProteinDetectionProtocolType_Properties_SetAndGetA()
        {
            var obj = new ProteinDetectionProtocolType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC9" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type9" } };

            obj.id = "pdp1";
            obj.analysisSoftware_ref = "as2";

            Assert.That(obj.id, Is.EqualTo("pdp1"));
            Assert.That(obj.analysisSoftware_ref, Is.EqualTo("as2"));
        }

        [Test]
        public void SearchDatabaseType_Properties_SetAndGet()
        {
            var obj = new SearchDatabaseType();
            var cvParams = new CVParamType[] { new CVParamType { accession = "ACC1" } };
            var userParams = new UserParamType[] { new UserParamType { type = "type1" } };

            obj.id = "db1";
            obj.location = "file://db.fasta";
            obj.name = "TestDB";
            obj.numDatabaseSequences = 1234;
            obj.numResidues = 56789;
            obj.cvParam = cvParams;

            Assert.That(obj.id, Is.EqualTo("db1"));
            Assert.That(obj.location, Is.EqualTo("file://db.fasta"));
            Assert.That(obj.name, Is.EqualTo("TestDB"));
            Assert.That(obj.numDatabaseSequences, Is.EqualTo(1234));
            Assert.That(obj.numResidues, Is.EqualTo(56789));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
        }

        [Test]
        public void InputsType_Properties_SetAndGet()
        {
            var obj = new InputsType();
            var sourceFiles = new SourceFileType[] { new SourceFileType() };
            var searchDatabases = new SearchDatabaseType[] { new SearchDatabaseType() };
            var spectraData = new SpectraDataType[] { new SpectraDataType() };

            obj.SourceFile = sourceFiles;
            obj.SearchDatabase = searchDatabases;
            obj.SpectraData = spectraData;

            Assert.That(obj.SourceFile, Is.EqualTo(sourceFiles));
            Assert.That(obj.SearchDatabase, Is.EqualTo(searchDatabases));
            Assert.That(obj.SpectraData, Is.EqualTo(spectraData));
        }

        [Test]
        public void SourceFileType_Properties_SetAndGetB()
        {
            var obj = new SourceFileType();
            obj.id = "sf1";
            obj.location = "file://test.mzML";
            obj.name = "TestFile";

            Assert.That(obj.id, Is.EqualTo("sf1"));
            Assert.That(obj.location, Is.EqualTo("file://test.mzML"));
            Assert.That(obj.name, Is.EqualTo("TestFile"));
        }

        [Test]
        public void SearchDatabaseType_Properties_SetAndGetA()
        {
            var obj = new SearchDatabaseType();
            obj.id = "db1";
            obj.location = "file://db.fasta";
            obj.name = "TestDB";

            Assert.That(obj.id, Is.EqualTo("db1"));
            Assert.That(obj.location, Is.EqualTo("file://db.fasta"));
            Assert.That(obj.name, Is.EqualTo("TestDB"));
        }

        [Test]
        public void SpectraDataType_Properties_SetAndGetB()
        {
            var obj = new SpectraDataType();
            obj.id = "sd1";
            obj.location = "file://spectra.mzML";
            obj.name = "Spectra1";

            Assert.That(obj.id, Is.EqualTo("sd1"));
            Assert.That(obj.location, Is.EqualTo("file://spectra.mzML"));
            Assert.That(obj.name, Is.EqualTo("Spectra1"));
        }
        [Test]
        public void MzIdentMLType111_Properties_SetAndGet()
        {
            var obj = new MzIdentMLType111();
            var cvList = new cvType[] { new cvType { id = "cv1" } };
            var analysisSoftwareList = new AnalysisSoftwareType[] { new AnalysisSoftwareType() };
            var provider = new ProviderType();
            var auditCollection = new AbstractContactType[] { new PersonType() };
            var analysisSampleCollection = new SampleType[] { new SampleType() };
            var sequenceCollection = new SequenceCollectionType();
            var analysisCollection = new AnalysisCollectionType();
            var analysisProtocolCollection = new AnalysisProtocolCollectionType();
            var dataCollection = new DataCollectionType();
            var bibliographicReference = new BibliographicReferenceType[] { new BibliographicReferenceType() };
            var creationDate = new System.DateTime(2024, 1, 1);
            var version = "1.1.1";

            obj.cvList = cvList;
            obj.AnalysisSoftwareList = analysisSoftwareList;
            obj.Provider = provider;
            obj.AuditCollection = auditCollection;
            obj.AnalysisSampleCollection = analysisSampleCollection;
            obj.SequenceCollection = sequenceCollection;
            obj.AnalysisCollection = analysisCollection;
            obj.AnalysisProtocolCollection = analysisProtocolCollection;
            obj.DataCollection = dataCollection;
            obj.BibliographicReference = bibliographicReference;
            obj.creationDate = creationDate;
            obj.creationDateSpecified = true;
            obj.version = version;

            Assert.That(obj.cvList, Is.EqualTo(cvList));
            Assert.That(obj.AnalysisSoftwareList, Is.EqualTo(analysisSoftwareList));
            Assert.That(obj.Provider, Is.EqualTo(provider));
            Assert.That(obj.AuditCollection, Is.EqualTo(auditCollection));
            Assert.That(obj.AnalysisSampleCollection, Is.EqualTo(analysisSampleCollection));
            Assert.That(obj.SequenceCollection, Is.EqualTo(sequenceCollection));
            Assert.That(obj.AnalysisCollection, Is.EqualTo(analysisCollection));
            Assert.That(obj.AnalysisProtocolCollection, Is.EqualTo(analysisProtocolCollection));
            Assert.That(obj.DataCollection, Is.EqualTo(dataCollection));
            Assert.That(obj.BibliographicReference, Is.EqualTo(bibliographicReference));
            Assert.That(obj.creationDate, Is.EqualTo(creationDate));
            Assert.That(obj.creationDateSpecified, Is.True);
            Assert.That(obj.version, Is.EqualTo(version));
        }

        [Test]
        public void cvType_Properties_SetAndGetC()
        {
            var obj = new cvType();
            obj.fullName = "Full Name";
            obj.version = "1.0";
            obj.uri = "http://example.com";
            obj.id = "cv1";

            Assert.That(obj.fullName, Is.EqualTo("Full Name"));
            Assert.That(obj.version, Is.EqualTo("1.0"));
            Assert.That(obj.uri, Is.EqualTo("http://example.com"));
            Assert.That(obj.id, Is.EqualTo("cv1"));
        }

        [Test]
        public void SpectrumIdentificationItemRefType_Properties_SetAndGetC()
        {
            var obj = new SpectrumIdentificationItemRefType();
            obj.spectrumIdentificationItem_ref = "item1";
            Assert.That(obj.spectrumIdentificationItem_ref, Is.EqualTo("item1"));
        }

        [Test]
        public void PeptideHypothesisType_Properties_SetAndGetC()
        {
            var obj = new PeptideHypothesisType();
            var refs = new SpectrumIdentificationItemRefType[] { new SpectrumIdentificationItemRefType { spectrumIdentificationItem_ref = "item1" } };
            obj.SpectrumIdentificationItemRef = refs;
            obj.peptideEvidence_ref = "pepEv1";
            Assert.That(obj.SpectrumIdentificationItemRef, Is.EqualTo(refs));
            Assert.That(obj.peptideEvidence_ref, Is.EqualTo("pepEv1"));
        }

        [Test]
        public void FragmentArrayType_Properties_SetAndGet()
        {
            var obj = new FragmentArrayType();
            var values = new float[] { 1.1f, 2.2f };
            obj.values = values;
            obj.measure_ref = "m1";
            Assert.That(obj.values, Is.EqualTo(values));
            Assert.That(obj.measure_ref, Is.EqualTo("m1"));
        }

        [Test]
        public void IonTypeType_Properties_SetAndGet()
        {
            var obj = new IonTypeType();
            var fragmentArrays = new FragmentArrayType[] { new FragmentArrayType() };
            var cvParam = new CVParamType();
            var index = new string[] { "1", "2" };
            obj.FragmentArray = fragmentArrays;
            obj.cvParam = cvParam;
            obj.index = index;
            obj.charge = 2;
            Assert.That(obj.FragmentArray, Is.EqualTo(fragmentArrays));
            Assert.That(obj.cvParam, Is.EqualTo(cvParam));
            Assert.That(obj.index, Is.EqualTo(index));
            Assert.That(obj.charge, Is.EqualTo(2));
        }

        [Test]
        public void CVParamType_Properties_SetAndGetB()
        {
            var obj = new CVParamType();
            obj.cvRef = "cvRef1";
            obj.accession = "ACC1";
            Assert.That(obj.cvRef, Is.EqualTo("cvRef1"));
            Assert.That(obj.accession, Is.EqualTo("ACC1"));
        }


        [Test]
        public void UserParamType_Properties_SetAndGetB()
        {
            var obj = new UserParamType();
            obj.type = "type1";
            Assert.That(obj.type, Is.EqualTo("type1"));
        }

        [Test]
        public void PeptideEvidenceRefType_Properties_SetAndGetB()
        {
            var obj = new PeptideEvidenceRefType();
            obj.peptideEvidence_ref = "pepEv1";
            Assert.That(obj.peptideEvidence_ref, Is.EqualTo("pepEv1"));
        }

        [Test]
        public void AnalysisDataType_Properties_SetAndGetB()
        {
            var obj = new AnalysisDataType();
            var sil = new SpectrumIdentificationListType[] { new SpectrumIdentificationListType() };
            var pdl = new ProteinDetectionListType();
            obj.SpectrumIdentificationList = sil;
            obj.ProteinDetectionList = pdl;
            Assert.That(obj.SpectrumIdentificationList, Is.EqualTo(sil));
            Assert.That(obj.ProteinDetectionList, Is.EqualTo(pdl));
        }

        [Test]
        public void SpectrumIdentificationListType_Properties_SetAndGetB()
        {
            var obj = new SpectrumIdentificationListType();
            var fragTable = new MeasureType[] { new MeasureType() };
            var sir = new SpectrumIdentificationResultType[] { new SpectrumIdentificationResultType() };
            var cvParams = new CVParamType[] { new CVParamType() };
            var userParams = new UserParamType[] { new UserParamType() };
            obj.FragmentationTable = fragTable;
            obj.SpectrumIdentificationResult = sir;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.numSequencesSearched = 123;
            obj.numSequencesSearchedSpecified = true;
            Assert.That(obj.FragmentationTable, Is.EqualTo(fragTable));
            Assert.That(obj.SpectrumIdentificationResult, Is.EqualTo(sir));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.numSequencesSearched, Is.EqualTo(123));
            Assert.That(obj.numSequencesSearchedSpecified, Is.True);
        }

        [Test]
        public void MeasureType_Properties_SetAndGet()
        {
            var obj = new MeasureType();
            var cvParams = new CVParamType[] { new CVParamType() };
            obj.cvParam = cvParams;
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
        }
        [Test]
        public void SpectrumIdentificationResultType_Properties_SetAndGetB()
        {
            var obj = new SpectrumIdentificationResultType();
            var sii = new SpectrumIdentificationItemType[] { new SpectrumIdentificationItemType() };
            var cvParams = new CVParamType[] { new CVParamType() };
            var userParams = new UserParamType[] { new UserParamType() };

            obj.spectrumID = "sid1";
            obj.spectraData_ref = "sd1";
            obj.SpectrumIdentificationItem = sii;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.spectrumID, Is.EqualTo("sid1"));
            Assert.That(obj.spectraData_ref, Is.EqualTo("sd1"));
            Assert.That(obj.SpectrumIdentificationItem, Is.EqualTo(sii));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void SpectrumIdentificationItemType_Properties_SetAndGetB()
        {
            var obj = new SpectrumIdentificationItemType();
            var peptideEvidenceRefs = new PeptideEvidenceRefType[] { new PeptideEvidenceRefType() };
            var fragmentations = new IonTypeType[] { new IonTypeType() };
            var cvParams = new CVParamType[] { new CVParamType() };
            var userParams = new UserParamType[] { new UserParamType() };

            obj.id = "sii1";
            obj.chargeState = 2;
            obj.experimentalMassToCharge = 500.2;
            obj.calculatedMassToCharge = 499.8;
            obj.peptide_ref = "pep1";
            obj.rank = 1;
            obj.passThreshold = true;
            obj.PeptideEvidenceRef = peptideEvidenceRefs;
            obj.Fragmentation = fragmentations;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("sii1"));
            Assert.That(obj.chargeState, Is.EqualTo(2));
            Assert.That(obj.experimentalMassToCharge, Is.EqualTo(500.2));
            Assert.That(obj.calculatedMassToCharge, Is.EqualTo(499.8));
            Assert.That(obj.peptide_ref, Is.EqualTo("pep1"));
            Assert.That(obj.rank, Is.EqualTo(1));
            Assert.That(obj.passThreshold, Is.True);
            Assert.That(obj.PeptideEvidenceRef, Is.EqualTo(peptideEvidenceRefs));
            Assert.That(obj.Fragmentation, Is.EqualTo(fragmentations));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void ProteinDetectionListType_Properties_SetAndGetB()
        {
            var obj = new ProteinDetectionListType();
            var proteinAmbiguityGroups = new ProteinAmbiguityGroupType[] { new ProteinAmbiguityGroupType() };
            var cvParams = new CVParamType[] { new CVParamType() };
            var userParams = new UserParamType[] { new UserParamType() };

            obj.id = "pdl1";
            obj.ProteinAmbiguityGroup = proteinAmbiguityGroups;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("pdl1"));
            Assert.That(obj.ProteinAmbiguityGroup, Is.EqualTo(proteinAmbiguityGroups));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void ProteinAmbiguityGroupType_Properties_SetAndGetB()
        {
            var obj = new ProteinAmbiguityGroupType();
            var proteinDetectionHypotheses = new ProteinDetectionHypothesisType[] { new ProteinDetectionHypothesisType() };
            var cvParams = new CVParamType[] { new CVParamType() };
            var userParams = new UserParamType[] { new UserParamType() };

            obj.id = "pag1";
            obj.ProteinDetectionHypothesis = proteinDetectionHypotheses;
            obj.cvParam = cvParams;
            obj.userParam = userParams;

            Assert.That(obj.id, Is.EqualTo("pag1"));
            Assert.That(obj.ProteinDetectionHypothesis, Is.EqualTo(proteinDetectionHypotheses));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
        }

        [Test]
        public void ProteinDetectionHypothesisType_Properties_SetAndGetB()
        {
            var obj = new ProteinDetectionHypothesisType();
            var peptideHypotheses = new PeptideHypothesisType[] { new PeptideHypothesisType() };
            var cvParams = new CVParamType[] { new CVParamType() };
            var userParams = new UserParamType[] { new UserParamType() };

            obj.id = "pdh1";
            obj.dBSequence_ref = "dbs1";
            obj.PeptideHypothesis = peptideHypotheses;
            obj.cvParam = cvParams;
            obj.userParam = userParams;
            obj.passThreshold = true;

            Assert.That(obj.id, Is.EqualTo("pdh1"));
            Assert.That(obj.dBSequence_ref, Is.EqualTo("dbs1"));
            Assert.That(obj.PeptideHypothesis, Is.EqualTo(peptideHypotheses));
            Assert.That(obj.cvParam, Is.EqualTo(cvParams));
            Assert.That(obj.userParam, Is.EqualTo(userParams));
            Assert.That(obj.passThreshold, Is.True);
        }
        private string TestFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratible_Yeast.mzID");

        private MzIdentMLType111 LoadMzIdentML()
        {
            var serializer = new XmlSerializer(typeof(MzIdentMLType111));
            using (var stream = File.OpenRead(TestFilePath))
            {
                return (MzIdentMLType111)serializer.Deserialize(stream);
            }
        }

        [Test]
        public void CanDeserializeMzIdentMLFile()
        {
            var mzIdentML = LoadMzIdentML();
            Assert.That(mzIdentML, Is.Not.Null);
            Assert.That(string.IsNullOrEmpty(mzIdentML.version), Is.False);
            Assert.That(mzIdentML.creationDate > DateTime.MinValue, Is.True);
        }

        [Test]
        public void CvList_IsNotEmpty_And_HasExpectedFields()
        {
            var mzIdentML = LoadMzIdentML();
            Assert.That(mzIdentML.cvList, Is.Not.Null);
            Assert.That(mzIdentML.cvList.Length, Is.GreaterThan(0));

            foreach (var cv in mzIdentML.cvList)
            {
                Assert.That(string.IsNullOrEmpty(cv.id), Is.False);
                Assert.That(string.IsNullOrEmpty(cv.fullName), Is.False);
                Assert.That(string.IsNullOrEmpty(cv.uri), Is.False);
            }
        }

        [Test]
        public void AnalysisSoftwareList_ContainsExpectedSoftware()
        {
            var mzIdentML = LoadMzIdentML();
            Assert.That(mzIdentML.AnalysisSoftwareList, Is.Not.Null);
            Assert.That(mzIdentML.AnalysisSoftwareList.Length, Is.GreaterThan(0));
        }

        [Test]
        public void DataCollection_And_AnalysisData_ArePresent()
        {
            var mzIdentML = LoadMzIdentML();
            Assert.That(mzIdentML.DataCollection, Is.Not.Null);
            Assert.That(mzIdentML.DataCollection.AnalysisData, Is.Not.Null);
            Assert.That(mzIdentML.DataCollection.AnalysisData.SpectrumIdentificationList, Is.Not.Null);
            Assert.That(mzIdentML.DataCollection.AnalysisData.SpectrumIdentificationList.Length, Is.GreaterThan(0));
        }

        [Test]
        public void SpectrumIdentificationResult_HasItems()
        {
            var mzIdentML = LoadMzIdentML();
            var sil = mzIdentML.DataCollection.AnalysisData.SpectrumIdentificationList[0];
            Assert.That(sil.SpectrumIdentificationResult, Is.Not.Null);
            Assert.That(sil.SpectrumIdentificationResult.Length, Is.GreaterThan(0));

            foreach (var result in sil.SpectrumIdentificationResult)
            {
                Assert.That(string.IsNullOrEmpty(result.spectrumID), Is.False);
                Assert.That(result.SpectrumIdentificationItem, Is.Not.Null);
                Assert.That(result.SpectrumIdentificationItem.Length, Is.GreaterThan(0));
            }
        }

        [Test]
        public void ProteinDetectionList_And_ProteinAmbiguityGroup_ArePresent()
        {
            var mzIdentML = LoadMzIdentML();
            var pdl = mzIdentML.DataCollection.AnalysisData.ProteinDetectionList;
            Assert.That(pdl, Is.Not.Null);
            Assert.That(pdl.ProteinAmbiguityGroup, Is.Not.Null);
            Assert.That(pdl.ProteinAmbiguityGroup.Length, Is.GreaterThan(0));

            foreach (var pag in pdl.ProteinAmbiguityGroup)
            {
                Assert.That(pag.ProteinDetectionHypothesis, Is.Not.Null);
                Assert.That(pag.ProteinDetectionHypothesis.Length, Is.GreaterThan(0));
            }
        }


        [Test]
        public void CvList_ContainsExpectedCvRefsAndAccessions()
        {
            var mzIdentML = LoadMzIdentML();
            Assert.That(mzIdentML.cvList, Is.Not.Null);
            Assert.That(mzIdentML.cvList.Length, Is.EqualTo(4));
            Assert.That(mzIdentML.cvList[0].id, Is.EqualTo("PSI-MS"));
            Assert.That(mzIdentML.cvList[0].fullName, Is.EqualTo("Proteomics Standards Initiative Mass Spectrometry Vocabularies"));
            Assert.That(mzIdentML.cvList[1].id, Is.EqualTo("PSI-MOD"));
            Assert.That(mzIdentML.cvList[2].id, Is.EqualTo("UNIMOD"));
            Assert.That(mzIdentML.cvList[3].id, Is.EqualTo("UO"));
        }

        [Test]
        public void AnalysisSoftwareList_ContainsExpectedSoftwareC()
        {
            var mzIdentML = LoadMzIdentML();
            Assert.That(mzIdentML.AnalysisSoftwareList, Is.Not.Null);
            Assert.That(mzIdentML.AnalysisSoftwareList.Length, Is.EqualTo(1));
            Assert.That(mzIdentML.AnalysisSoftwareList[0].id, Is.EqualTo("AS_MetaMorpheus"));
            Assert.That(mzIdentML.AnalysisSoftwareList[0].name, Is.EqualTo("MetaMorpheus"));
            Assert.That(mzIdentML.AnalysisSoftwareList[0].version, Is.EqualTo("1.10.0"));
        }

        [Test]
        public void CvParam_ContainsExpectedCvRefAndAccession()
        {
            var mzIdentML = LoadMzIdentML();
            var software = mzIdentML.AnalysisSoftwareList[0];
            var roleCvParam = software.ContactRole.Role.cvParam;
            Assert.That(roleCvParam.cvRef, Is.EqualTo("PSI-MS"));
            Assert.That(roleCvParam.accession, Is.EqualTo("MS:1001267"));
        }
    }
}
