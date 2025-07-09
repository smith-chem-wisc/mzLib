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
        public void Property_GetSet_Works()
        {
            var obj = new MzIdentMLType111();

            var cvList = new cvType[] { new cvType { id = "cv1" } };
            obj.cvList = cvList;
            Assert.That(obj.cvList, Is.EqualTo(cvList));

            var analysisSoftwareList = new AnalysisSoftwareType[] { new AnalysisSoftwareType() };
            obj.AnalysisSoftwareList = analysisSoftwareList;
            Assert.That(obj.AnalysisSoftwareList, Is.EqualTo(analysisSoftwareList));

            var provider = new ProviderType();
            obj.Provider = provider;
            Assert.That(obj.Provider, Is.EqualTo(provider));

            var auditCollection = new AbstractContactType[] { new OrganizationType() };
            obj.AuditCollection = auditCollection;
            Assert.That(obj.AuditCollection, Is.EqualTo(auditCollection));

            var sampleCollection = new SampleType[] { new SampleType() };
            obj.AnalysisSampleCollection = sampleCollection;
            Assert.That(obj.AnalysisSampleCollection, Is.EqualTo(sampleCollection));

            var sequenceCollection = new SequenceCollectionType();
            obj.SequenceCollection = sequenceCollection;
            Assert.That(obj.SequenceCollection, Is.EqualTo(sequenceCollection));

            var analysisCollection = new AnalysisCollectionType();
            obj.AnalysisCollection = analysisCollection;
            Assert.That(obj.AnalysisCollection, Is.EqualTo(analysisCollection));

            var analysisProtocolCollection = new AnalysisProtocolCollectionType();
            obj.AnalysisProtocolCollection = analysisProtocolCollection;
            Assert.That(obj.AnalysisProtocolCollection, Is.EqualTo(analysisProtocolCollection));

            var dataCollection = new DataCollectionType();
            obj.DataCollection = dataCollection;
            Assert.That(obj.DataCollection, Is.EqualTo(dataCollection));

            var biblioRefs = new BibliographicReferenceType[] { new BibliographicReferenceType() };
            obj.BibliographicReference = biblioRefs;
            Assert.That(obj.BibliographicReference, Is.EqualTo(biblioRefs));

            var date = new DateTime(2020, 1, 1);
            obj.creationDate = date;
            Assert.That(obj.creationDate, Is.EqualTo(date));

            obj.creationDateSpecified = true;
            Assert.That(obj.creationDateSpecified, Is.EqualTo(true));

            obj.version = "1.1.1";
            Assert.That(obj.version, Is.EqualTo("1.1.1"));
        }
    }

    [TestFixture]
    public class cvTypeTests
    {
        [Test]
        public void Property_GetSet_Works()
        {
            var obj = new cvType();

            obj.fullName = "Full Name";
            Assert.That(obj.fullName, Is.EqualTo("Full Name"));

            obj.version = "1.0";
            Assert.That(obj.version, Is.EqualTo("1.0"));

            obj.uri = "http://example.com";
            Assert.That(obj.uri, Is.EqualTo("http://example.com"));

            obj.id = "cv1";
            Assert.That(obj.id, Is.EqualTo("cv1"));
        }
    }

    [TestFixture]
    public class SpectrumIdentificationItemRefTypeTests
    {
        [Test]
        public void Property_GetSet_Works()
        {
            var obj = new SpectrumIdentificationItemRefType();
            obj.spectrumIdentificationItem_ref = "ref1";
            Assert.That(obj.spectrumIdentificationItem_ref, Is.EqualTo("ref1"));
        }
    }

    [TestFixture]
    public class PeptideHypothesisTypeTests
    {
        [Test]
        public void Property_GetSet_Works()
        {
            var obj = new PeptideHypothesisType();

            var refs = new SpectrumIdentificationItemRefType[] { new SpectrumIdentificationItemRefType() };
            obj.SpectrumIdentificationItemRef = refs;
            Assert.That(obj.SpectrumIdentificationItemRef, Is.EqualTo(refs));

            obj.peptideEvidence_ref = "pepRef";
            Assert.That(obj.peptideEvidence_ref, Is.EqualTo("pepRef"));
        }
    }

    [TestFixture]
    public class FragmentArrayTypeTests
    {
        [Test]
        public void Property_GetSet_Works()
        {
            var obj = new FragmentArrayType();

            var values = new float[] { 1.1f, 2.2f };
            obj.values = values;
            Assert.That(obj.values, Is.EqualTo(values));

            obj.measure_ref = "measure1";
            Assert.That(obj.measure_ref, Is.EqualTo("measure1"));
        }
    }

    [TestFixture]
    public class IonTypeTypeTests
    {
        [Test]
        public void Property_GetSet_Works()
        {
            var obj = new IonTypeType();

            var fragmentArray = new FragmentArrayType[] { new FragmentArrayType() };
            obj.FragmentArray = fragmentArray;
            Assert.That(obj.FragmentArray, Is.EqualTo(fragmentArray));

            var cvParam = new CVParamType();
            obj.cvParam = cvParam;
            Assert.That(obj.cvParam, Is.EqualTo(cvParam));

            var index = new string[] { "1", "2" };
            obj.index = index;
            Assert.That(obj.index, Is.EqualTo(index));

            obj.charge = 2;
            Assert.That(obj.charge, Is.EqualTo(2));
        }
    }

    [TestFixture]
    public class CVParamTypeTests
    {
        [Test]
        public void Property_GetSet_Works()
        {
            var obj = new CVParamType();

            obj.cvRef = "cvRef";
            Assert.That(obj.cvRef, Is.EqualTo("cvRef"));

            obj.accession = "acc";
            Assert.That(obj.accession, Is.EqualTo("acc"));
        }
    }
}
