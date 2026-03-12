using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using System.Collections.Generic;
using System.Linq;
using static Test.Omics.SequenceConversion.GroundTruthTestData;

namespace Test.Omics.SequenceConversion;

/// <summary>
/// Tests for serializing canonical sequences into various format strings.
/// </summary>
[TestFixture]
public class SerializationTests
{
    private MzLibSequenceSerializer _mzLibSerializer;
    private MassShiftSequenceSerializer _massShiftSerializer;
    private ChronologerSequenceSerializer _chronologerSerializer;
    private UnimodSequenceSerializer _unimodSerializer;
    private UniProtSequenceSerializer _uniprotSerializer;
    private UniProtSequenceSerializer _uniprotSerializerWithModType;
    private MzLibSequenceParser _mzLibParser;
    private MassShiftSequenceParser _massShiftParser;

    public static IEnumerable<UniProtMappingTestCase> UniProtMappingTestCases() => GroundTruthTestData.UniProtMappingTestCases;
    public static IEnumerable<SequenceConversionTestCase> CoreTestCases() => GroundTruthTestData.CoreTestCases;
    public static IEnumerable<SequenceConversionTestCase> EdgeCases() => GroundTruthTestData.EdgeCases;

    [SetUp]
    public void Setup()
    {
        _mzLibSerializer = new MzLibSequenceSerializer();
        _massShiftSerializer = new MassShiftSequenceSerializer(new(4));
        _chronologerSerializer = new ChronologerSequenceSerializer();
        _unimodSerializer = new UnimodSequenceSerializer();
        _uniprotSerializer = new UniProtSequenceSerializer();
        _uniprotSerializerWithModType = new UniProtSequenceSerializer(new UniProtSequenceSchema(true, false));
        _mzLibParser = new MzLibSequenceParser();
        _massShiftParser = new MassShiftSequenceParser();
    }

    [Test]
    [TestCaseSource(nameof(CoreTestCases))]
    public void MzLibSerializer_CoreTestCases_SerializesCorrectly(SequenceConversionTestCase testCase)
    {
        // Arrange - parse to get canonical form
        var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
        Assert.That(canonical, Is.Not.Null);

        // Act
        var result = _mzLibSerializer.Serialize(canonical.Value);

        // Assert
        Assert.That(result, Is.Not.Null);
        Assert.That(result, Is.EqualTo(testCase.MzLibFormat));
    }

    [Test]
    [TestCaseSource(nameof(CoreTestCases))]
    public void ChronologerSerializer_CoreTestCases_SerializesCorrectly(SequenceConversionTestCase testCase)
    {
        // Arrange
        var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
        Assert.That(canonical, Is.Not.Null);

        // Act
        var result = _chronologerSerializer.Serialize(canonical.Value, null, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        // Assert
        Assert.That(result, Is.Not.Null);
        Assert.That(result, Is.EqualTo(testCase.ChronologerFormat));
    }

    [Test]
    [TestCaseSource(nameof(CoreTestCases))]
    public void UnimodSerializer_CoreTestCases_ConvertsCorrectly(SequenceConversionTestCase testCase)
    {
        // Arrange - parse from mzLib format
        var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
        Assert.That(canonical, Is.Not.Null);

        // Act - serialize to Unimod format
        var result = _unimodSerializer.Serialize(canonical.Value);

        // Assert
        Assert.That(result, Is.Not.Null);
        Assert.That(result, Is.EqualTo(testCase.UnimodUpperCaseFormat));
    }

    [Test]
    [TestCaseSource(nameof(CoreTestCases))]
    public void UniprotSerializer_KnownCases_SerializesCorrectly(SequenceConversionTestCase testCase)
    {
        var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
        Assert.That(canonical, Is.Not.Null);

        var result = _uniprotSerializer.Serialize(canonical.Value);

        Assert.That(result, Is.Not.Null);
        Assert.That(result, Is.EqualTo(testCase.UniProtFormat));
    }

    /// <summary>
    /// Tests that specific mzLib modifications correctly map to their UniProt representations.
    /// Each test case validates one modification type on a specific residue.
    /// </summary>
    [Test]
    [TestCaseSource(nameof(UniProtMappingTestCases))]
    public void UniprotSerializer_ModificationMappings_ResolveCorrectly(UniProtMappingTestCase testCase)
    {
        // Build a simple test sequence with the modification
        string mzLibSequence;
        string expectedUniProtSequence;

        if (testCase.IsNTerminal)
        {
            // N-terminal modification: [mod]XEPTIDE where X is the target residue
            mzLibSequence = $"[{testCase.MzLibModification}]{testCase.TargetResidue}EPTIDE";
            expectedUniProtSequence = $"[{testCase.ExpectedUniProtName}]{testCase.TargetResidue}EPTIDE";
        }
        else
        {
            // Internal modification: PEP{X}[mod]IDE where X is the target residue
            mzLibSequence = $"PEP{testCase.TargetResidue}[{testCase.MzLibModification}]IDE";
            expectedUniProtSequence = $"PEP{testCase.TargetResidue}[{testCase.ExpectedUniProtName}]IDE";
        }

        // Parse the mzLib sequence
        var canonical = _mzLibParser.Parse(mzLibSequence);
        Assert.That(canonical, Is.Not.Null, $"Failed to parse mzLib sequence: {mzLibSequence}");

        // Serialize to UniProt format
        var result = _uniprotSerializer.Serialize(canonical.Value);

        // Verify the result
        Assert.That(result, Is.Not.Null, $"UniProt serialization returned null for: {mzLibSequence}");
        Assert.That(result, Is.EqualTo(expectedUniProtSequence), 
            $"Expected UniProt format '{expectedUniProtSequence}' but got '{result}'");
    }


    [Test]
    [TestCaseSource(nameof(UniProtMappingTestCases))]
    public void UniprotSerializerWithModType_ModificationMappings_ResolveCorrectly(UniProtMappingTestCase testCase)
    {
        // Build a simple test sequence with the modification
        string mzLibSequence;
        string expectedUniProtSequence;

        if (testCase.IsNTerminal)
        {
            // N-terminal modification: [mod]XEPTIDE where X is the target residue
            mzLibSequence = $"[{testCase.MzLibModification}]{testCase.TargetResidue}EPTIDE";
            expectedUniProtSequence = $"[UniProt:{testCase.ExpectedUniProtName}]{testCase.TargetResidue}EPTIDE";
        }
        else
        {
            // Internal modification: PEP{X}[mod]IDE where X is the target residue
            mzLibSequence = $"PEP{testCase.TargetResidue}[{testCase.MzLibModification}]IDE";
            expectedUniProtSequence = $"PEP{testCase.TargetResidue}[UniProt:{testCase.ExpectedUniProtName}]IDE";
        }

        // Parse the mzLib sequence
        var canonical = _mzLibParser.Parse(mzLibSequence);
        Assert.That(canonical, Is.Not.Null, $"Failed to parse mzLib sequence: {mzLibSequence}");

        // Serialize to UniProt format
        var result = _uniprotSerializerWithModType.Serialize(canonical.Value);

        // Verify the result
        Assert.That(result, Is.Not.Null, $"UniProt serialization returned null for: {mzLibSequence}");
        Assert.That(result, Is.EqualTo(expectedUniProtSequence),
            $"Expected UniProt format '{expectedUniProtSequence}' but got '{result}'");
    }

    [Test]
    public void MzLibSerializer_EmptySequence_ReturnsNull()
    {
        // Arrange
        var canonical = CanonicalSequence.Empty;

        // Act
        var result = _mzLibSerializer.Serialize(canonical, null, SequenceConversionHandlingMode.ReturnNull);

        // Assert
        Assert.That(result, Is.Null);
    }

    #region MassShift Serializer Tests

    [Test]
    [TestCaseSource(nameof(CoreTestCases))]
    public void MassShiftSerializer_CoreTestCases_SerializesCorrectly(SequenceConversionTestCase testCase)
    {
        // Arrange - parse from MassShift format to get canonical form
        var canonical = _massShiftParser.Parse(testCase.MassShiftFormat);
        Assert.That(canonical, Is.Not.Null);

        // Act - serialize back to MassShift format
        var result = _massShiftSerializer.Serialize(canonical.Value);

        // Assert
        Assert.That(result, Is.Not.Null);
        Assert.That(result, Is.EqualTo(testCase.MassShiftFormat));
    }

    [Test]
    [TestCaseSource(nameof(EdgeCases))]
    public void MassShiftSerializer_EdgeCases_SerializesCorrectly(SequenceConversionTestCase testCase)
    {
        // Arrange
        var canonical = _massShiftParser.Parse(testCase.MassShiftFormat);
        Assert.That(canonical, Is.Not.Null);

        // Act
        var result = _massShiftSerializer.Serialize(canonical.Value);

        // Assert
        Assert.That(result, Is.Not.Null);
        Assert.That(result, Is.EqualTo(testCase.MassShiftFormat));
    }

    [Test]
    [TestCaseSource(nameof(CoreTestCases))]
    public void MzLibToMassShift_CoreTestCases_ConvertsCorrectly(SequenceConversionTestCase testCase)
    {
        // Arrange - parse from mzLib format
        var canonical = _mzLibParser.Parse(testCase.MzLibFormat);
        Assert.That(canonical, Is.Not.Null);

        // Act - serialize to MassShift format
        var result = _massShiftSerializer.Serialize(canonical.Value);

        // Assert
        Assert.That(result, Is.Not.Null);
        Assert.That(result, Is.EqualTo(testCase.MassShiftFormat));
    }

    [Test]
    public void MassShiftSerializer_EmptySequence_ReturnsNull()
    {
        // Arrange
        var canonical = CanonicalSequence.Empty;

        // Act
        var result = _massShiftSerializer.Serialize(canonical, null, SequenceConversionHandlingMode.ReturnNull);

        // Assert
        Assert.That(result, Is.Null);
    }

    [Test]
    public void MassShiftToMzLib_UsesStrictTypeAndIdWithMotifToken()
    {
        // Arrange
        var canonical = _massShiftParser.Parse("PEPTM[+15.9949]IDE");
        Assert.That(canonical, Is.Not.Null);

        // Act
        var result = _mzLibSerializer.Serialize(canonical.Value);

        // Assert
        Assert.That(result, Is.Not.Null);

        var openBracket = result!.IndexOf('[');
        var closeBracket = result.IndexOf(']', openBracket + 1);
        Assert.That(openBracket, Is.GreaterThanOrEqualTo(0));
        Assert.That(closeBracket, Is.GreaterThan(openBracket));

        var token = result.Substring(openBracket + 1, closeBracket - openBracket - 1);
        Assert.That(token, Does.Contain(":"));
        Assert.That(token.StartsWith(":"), Is.False);
        Assert.That(token.EndsWith(":"), Is.False);
    }

    [Test]
    [TestCase(UnimodLabelStyle.UpperCase, "UNIMOD")]
    [TestCase(UnimodLabelStyle.CamelCase, "Unimod")]
    [TestCase(UnimodLabelStyle.LowerCase, "unimod")]
    [TestCase(UnimodLabelStyle.NoLabel, "")]
    public void UnimodSerializer_LabelStyle_WritesExpectedToken(UnimodLabelStyle labelStyle, string expectedLabel)
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(2, "UNIMOD:35", unimodId: 35)
            .Build();

        var serializer = new UnimodSequenceSerializer(new UnimodSequenceFormatSchema(labelStyle));

        var result = serializer.Serialize(canonical);

        Assert.That(result, Is.Not.Null);
        var expectedToken = string.IsNullOrEmpty(expectedLabel) ? "[35]" : $"[{expectedLabel}:35]";
        Assert.That(result, Does.Contain(expectedToken));
    }

    [Test]
    public void UnimodSerializer_FromMzLib_ResolvesUnimodId()
    {
        var canonical = _mzLibParser.Parse("PEPTM[Common Variable:Oxidation on M]IDE");
        Assert.That(canonical, Is.Not.Null);

        var result = _unimodSerializer.Serialize(canonical.Value);

        Assert.That(result, Is.EqualTo("PEPTM[UNIMOD:35]IDE"));
    }

    [Test]
    public void EssentialSerializer_PrunesNonWhitelistedModificationTypes()
    {
        var canonical = _mzLibParser.Parse("[Common Biological:Acetylation on X]PEPTM[Common Variable:Oxidation on M]IDE");
        Assert.That(canonical, Is.Not.Null);

        var serializer = new EssentialSequenceSerializer(new Dictionary<string, int>
        {
            { "Common Variable", 0 }
        });

        var result = serializer.Serialize(canonical.Value);

        Assert.That(result, Is.EqualTo("PEPTM[Common Variable:Oxidation on M]IDE"));
    }

    #endregion
}
