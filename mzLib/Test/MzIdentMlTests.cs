using NUnit.Framework;
using System;
using System.Diagnostics.CodeAnalysis;
using mzIdentML111.Generated;

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
    }
}
