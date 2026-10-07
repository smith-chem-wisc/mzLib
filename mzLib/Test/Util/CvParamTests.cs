using MzLibUtil;
using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;

namespace Test.Util;

[TestFixture]
[ExcludeFromCodeCoverage]
public class CvParamTests
{
    [Test]
    public void ParameterlessConstructor_DefaultsAllPropertiesToEmptyString()
    {
        var cv = new CvParam();
        Assert.That(cv.CvLabel, Is.EqualTo(string.Empty));
        Assert.That(cv.Accession, Is.EqualTo(string.Empty));
        Assert.That(cv.Name, Is.EqualTo(string.Empty));
        Assert.That(cv.Value, Is.EqualTo(string.Empty));
        Assert.That(cv.UnitCvLabel, Is.EqualTo(string.Empty));
        Assert.That(cv.UnitAccession, Is.EqualTo(string.Empty));
        Assert.That(cv.UnitName, Is.EqualTo(string.Empty));
    }

    [Test]
    public void Constructor_SetsCoreFields_AndDefaultsUnitFieldsToEmpty()
    {
        var cv = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://ftp.pride.ebi.ac.uk/x");
        Assert.That(cv.CvLabel, Is.EqualTo("PRIDE"));
        Assert.That(cv.Accession, Is.EqualTo("PRIDE:0000469"));
        Assert.That(cv.Name, Is.EqualTo("FTP Protocol"));
        Assert.That(cv.Value, Is.EqualTo("ftp://ftp.pride.ebi.ac.uk/x"));
        // units are optional and absent here
        Assert.That(cv.UnitCvLabel, Is.EqualTo(string.Empty));
        Assert.That(cv.UnitAccession, Is.EqualTo(string.Empty));
        Assert.That(cv.UnitName, Is.EqualTo(string.Empty));
    }

    [Test]
    public void Constructor_SetsFullPsiFieldSetIncludingUnits()
    {
        var cv = new CvParam("MS", "MS:1000040", "m/z", "1234.56",
            unitCvLabel: "MS", unitAccession: "MS:1000040", unitName: "m/z");
        Assert.That(cv.UnitCvLabel, Is.EqualTo("MS"));
        Assert.That(cv.UnitAccession, Is.EqualTo("MS:1000040"));
        Assert.That(cv.UnitName, Is.EqualTo("m/z"));
    }

    [Test]
    public void ObjectInitializer_SetsAllProperties()
    {
        var cv = new CvParam
        {
            CvLabel = "PRIDE",
            Accession = "PRIDE:0000408",
            Name = "Search engine output file URI",
            Value = "SEARCH"
        };
        Assert.That(cv.Accession, Is.EqualTo("PRIDE:0000408"));
        Assert.That(cv.Value, Is.EqualTo("SEARCH"));
    }

    [Test]
    public void Equality_SameValues_AreEqual()
    {
        var a = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://x");
        var b = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://x");
        Assert.That(a, Is.EqualTo(b));
        Assert.That(a.GetHashCode(), Is.EqualTo(b.GetHashCode()));
    }

    [Test]
    public void Equality_DifferentAccession_AreNotEqual()
    {
        var a = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://x");
        var b = new CvParam("PRIDE", "PRIDE:0000468", "Aspera Protocol", "fasp://x");
        Assert.That(a, Is.Not.EqualTo(b));
    }

    [Test]
    public void EqualityOperators_CompareByValue()
    {
        var a = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://x");
        var b = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://x");
        var c = new CvParam("PRIDE", "PRIDE:0000468", "Aspera Protocol", "fasp://x");
        Assert.That(a == b, Is.True);
        Assert.That(a != c, Is.True);
    }

    [Test]
    public void With_ProducesModifiedCopy_LeavingOriginalUnchanged()
    {
        var a = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://x");
        var b = a with { Value = "ftp://y" };
        Assert.That(b.Accession, Is.EqualTo("PRIDE:0000469"));
        Assert.That(b.Value, Is.EqualTo("ftp://y"));
        Assert.That(a.Value, Is.EqualTo("ftp://x"));
    }

    [Test]
    public void ToString_IncludesAccessionAndName()
    {
        var cv = new CvParam("PRIDE", "PRIDE:0000469", "FTP Protocol", "ftp://x");
        var text = cv.ToString();
        Assert.That(text, Does.Contain("PRIDE:0000469"));
        Assert.That(text, Does.Contain("FTP Protocol"));
    }
}
