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
    private MzLibSequenceParser _mzLibParser;
    private MassShiftSequenceParser _massShiftParser;

    public static IEnumerable<SequenceConversionTestCase> CoreTestCases() => GroundTruthTestData.CoreTestCases;
    public static IEnumerable<SequenceConversionTestCase> EdgeCases() => GroundTruthTestData.EdgeCases;

    [SetUp]
    public void Setup()
    {
        _mzLibSerializer = new MzLibSequenceSerializer();
        _massShiftSerializer = new MassShiftSequenceSerializer(new(4));
        _chronologerSerializer = new ChronologerSequenceSerializer();
        _unimodSerializer = new UnimodSequenceSerializer();
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
    public void MassShiftSerializer_ModificationWithoutMass_RemoveModeSkipsModification()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(2, "Unknown:NoMass")
            .Build();
        var warnings = new ConversionWarnings();

        var result = _massShiftSerializer.Serialize(canonical, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);

        Assert.That(result, Is.EqualTo("PEPTIDE"));
        Assert.That(warnings.HasWarnings, Is.True);
        Assert.That(warnings.HasIncompatibleItems, Is.True);
    }

    [Test]
    public void MassShiftSerializer_ModificationWithoutMass_ThrowModeThrows()
    {
        var canonical = new CanonicalSequenceBuilder("PEPTIDE")
            .AddResidueModification(2, "Unknown:NoMass")
            .Build();

        Assert.That(
            () => _massShiftSerializer.Serialize(canonical, null, SequenceConversionHandlingMode.ThrowException),
            Throws.TypeOf<SequenceConversionException>());
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
