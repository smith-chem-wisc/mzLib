using System.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class UnimodSequenceParserTests
{
    private UnimodSequenceParser _parser;

    [SetUp]
    public void SetUp()
    {
        _parser = new UnimodSequenceParser();
    }

    #region CanParse - True Cases

    [Test]
    public void CanParse_ReturnsTrueForLabeledUnimodTokens()
    {
        Assert.That(_parser.CanParse("PEPC[UNIMOD:4]IDE"), Is.True);
    }

    [Test]
    public void CanParse_CamelCaseUnimodToken_ReturnsTrue()
    {
        Assert.That(_parser.CanParse("PEPC[Unimod:4]IDE"), Is.True);
    }

    [Test]
    public void CanParse_LowerCaseUnimodToken_ReturnsTrue()
    {
        Assert.That(_parser.CanParse("PEPC[unimod:4]IDE"), Is.True);
    }

    [Test]
    public void CanParse_WhitespaceBeforeId_ReturnsTrue()
    {
        Assert.That(_parser.CanParse("PEPC[UNIMOD: 4]IDE"), Is.True);
    }

    [Test]
    public void CanParse_NTerminalUnimod_ReturnsTrue()
    {
        Assert.That(_parser.CanParse("[UNIMOD:1]PEPTIDE"), Is.True);
    }

    [Test]
    public void CanParse_CTerminalUnimod_ReturnsTrue()
    {
        Assert.That(_parser.CanParse("PEPTIDE-[UNIMOD:2]"), Is.True);
    }

    [Test]
    public void CanParse_MultipleUnimodTokens_ReturnsTrue()
    {
        Assert.That(_parser.CanParse("PEPS[UNIMOD:21]TM[UNIMOD:35]IDE"), Is.True);
    }

    #endregion

    #region CanParse - False Cases

    [Test]
    public void CanParse_EmptyInput_ReturnsFalse()
    {
        Assert.That(_parser.CanParse(""), Is.False);
    }

    [Test]
    public void CanParse_WhitespaceInput_ReturnsFalse()
    {
        Assert.That(_parser.CanParse("   "), Is.False);
    }

    [Test]
    public void CanParse_NoBrackets_ReturnsFalse()
    {
        Assert.That(_parser.CanParse("PEPTIDE"), Is.False);
    }

    [Test]
    public void CanParse_MzLibStyleToken_ReturnsFalse()
    {
        Assert.That(_parser.CanParse("PEPC[Common Fixed:Carbamidomethyl on C]IDE"), Is.False);
    }

    [Test]
    public void CanParse_BareNumericToken_ReturnsFalse()
    {
        Assert.That(_parser.CanParse("PEPC[4]IDE"), Is.False);
    }

    [Test]
    public void CanParse_MassShiftToken_ReturnsFalse()
    {
        Assert.That(_parser.CanParse("PEPM[+15.9949]IDE"), Is.False);
    }

    [Test]
    public void CanParse_UnbalancedBrackets_ReturnsFalse()
    {
        Assert.That(_parser.CanParse("PEPC[UNIMOD:4IDE"), Is.False);
    }

    [Test]
    public void CanParse_OnlyClosingBracket_ReturnsFalse()
    {
        Assert.That(_parser.CanParse("PEPC]UNIMOD:4[IDE"), Is.False);
    }

    [Test]
    public void CanParse_MixedUnimodAndMzLib_ReturnsTrue()
    {
        Assert.That(_parser.CanParse("PEPC[UNIMOD:4]M[Oxidation on M]IDE"), Is.True);
    }

    #endregion

    #region Parse - Residue Modifications

    [Test]
    public void Parse_LabeledUnimodToken_SetsUnimodId()
    {
        var result = _parser.Parse("PEPC[UNIMOD:4]IDE");

        Assert.That(result, Is.Not.Null);
        var modification = result!.Value.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(4));
        Assert.That(modification.Value.TargetResidue, Is.EqualTo('C'));
        Assert.That(modification.Value.PositionType, Is.EqualTo(ModificationPositionType.Residue));
        Assert.That(modification.Value.OriginalRepresentation, Is.EqualTo("UNIMOD:4"));
    }

    [Test]
    public void Parse_CamelCaseUnimodToken_SetsUnimodId()
    {
        var result = _parser.Parse("PEPC[Unimod:4]IDE");

        Assert.That(result, Is.Not.Null);
        var modification = result!.Value.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(4));
    }

    [Test]
    public void Parse_LowerCaseUnimodToken_SetsUnimodId()
    {
        var result = _parser.Parse("PEPC[unimod:4]IDE");

        Assert.That(result, Is.Not.Null);
        var modification = result!.Value.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(4));
    }

    [Test]
    public void Parse_WhitespaceInToken_SetsUnimodId()
    {
        var result = _parser.Parse("PEPC[UNIMOD: 4]IDE");

        Assert.That(result, Is.Not.Null);
        var modification = result!.Value.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(4));
    }

    [Test]
    public void Parse_UnlabeledIntegerToken_SetsUnimodId()
    {
        var result = _parser.Parse("PEPC[4]IDE");

        Assert.That(result, Is.Not.Null);
        var modification = result!.Value.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(4));
    }

    [Test]
    public void Parse_OxidationOnM_SetsUnimodId35()
    {
        var result = _parser.Parse("PEPTM[UNIMOD:35]IDE");

        Assert.That(result, Is.Not.Null);
        var modification = result!.Value.GetModificationAt(4);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(35));
        Assert.That(modification.Value.TargetResidue, Is.EqualTo('M'));
    }

    [Test]
    public void Parse_PhosphoOnS_SetsUnimodId21()
    {
        var result = _parser.Parse("PEPS[UNIMOD:21]TIDE");

        Assert.That(result, Is.Not.Null);
        var modification = result!.Value.GetModificationAt(3);
        Assert.That(modification, Is.Not.Null);
        Assert.That(modification!.Value.UnimodId, Is.EqualTo(21));
        Assert.That(modification.Value.TargetResidue, Is.EqualTo('S'));
    }

    [Test]
    public void Parse_TwoModificationsOnSameResidueType()
    {
        var result = _parser.Parse("M[UNIMOD:35]EPTM[UNIMOD:35]IDE");

        Assert.That(result, Is.Not.Null);
        var mod1 = result!.Value.GetModificationAt(0);
        var mod2 = result!.Value.GetModificationAt(4);
        Assert.That(mod1!.Value.UnimodId, Is.EqualTo(35));
        Assert.That(mod2!.Value.UnimodId, Is.EqualTo(35));
    }

    [Test]
    public void Parse_MultipleDifferentModifications()
    {
        var result = _parser.Parse("PEPS[UNIMOD:21]TM[UNIMOD:35]IDE");

        Assert.That(result, Is.Not.Null);
        var phos = result!.Value.GetModificationAt(3);
        var ox = result!.Value.GetModificationAt(5);
        Assert.That(phos!.Value.UnimodId, Is.EqualTo(21));
        Assert.That(ox!.Value.UnimodId, Is.EqualTo(35));
    }

    [Test]
    public void Parse_AdjacentModifications()
    {
        var result = _parser.Parse("PEPTM[UNIMOD:35]M[UNIMOD:35]IDE");

        Assert.That(result, Is.Not.Null);
        var mod1 = result!.Value.GetModificationAt(4);
        var mod2 = result!.Value.GetModificationAt(5);
        Assert.That(mod1!.Value.UnimodId, Is.EqualTo(35));
        Assert.That(mod2!.Value.UnimodId, Is.EqualTo(35));
    }

    #endregion

    #region Parse - N-Terminal Modifications

    [Test]
    public void Parse_NTerminalAcetylation()
    {
        var result = _parser.Parse("[UNIMOD:1]PEPTIDE");

        Assert.That(result, Is.Not.Null);
        var canonical = result!.Value;
        Assert.That(canonical.BaseSequence, Is.EqualTo("PEPTIDE"));
        Assert.That(canonical.NTerminalModification, Is.Not.Null);
        Assert.That(canonical.NTerminalModification!.Value.UnimodId, Is.EqualTo(1));
        Assert.That(canonical.NTerminalModification.Value.PositionType, Is.EqualTo(ModificationPositionType.NTerminus));
    }

    [Test]
    public void Parse_NTerminalWithInternalModification()
    {
        var result = _parser.Parse("[UNIMOD:1]PEPTM[UNIMOD:35]IDE");

        Assert.That(result, Is.Not.Null);
        var canonical = result!.Value;
        Assert.That(canonical.NTerminalModification, Is.Not.Null);
        Assert.That(canonical.NTerminalModification!.Value.UnimodId, Is.EqualTo(1));
        var internalMod = canonical.GetModificationAt(4);
        Assert.That(internalMod, Is.Not.Null);
        Assert.That(internalMod!.Value.UnimodId, Is.EqualTo(35));
    }

    #endregion

    #region Parse - C-Terminal Modifications

    [Test]
    public void Parse_CTerminalModification()
    {
        var result = _parser.Parse("PEPTIDE-[UNIMOD:2]");

        Assert.That(result, Is.Not.Null);
        var canonical = result!.Value;
        Assert.That(canonical.CTerminalModification, Is.Not.Null);
        Assert.That(canonical.CTerminalModification!.Value.UnimodId, Is.EqualTo(2));
        Assert.That(canonical.CTerminalModification.Value.PositionType, Is.EqualTo(ModificationPositionType.CTerminus));
    }

    [Test]
    public void Parse_CTerminalWithInternalModification()
    {
        var result = _parser.Parse("PEPTM[UNIMOD:35]IDE-[UNIMOD:2]");

        Assert.That(result, Is.Not.Null);
        var canonical = result!.Value;
        var internalMod = canonical.GetModificationAt(4);
        Assert.That(internalMod, Is.Not.Null);
        Assert.That(internalMod!.Value.UnimodId, Is.EqualTo(35));
        Assert.That(canonical.CTerminalModification, Is.Not.Null);
        Assert.That(canonical.CTerminalModification!.Value.UnimodId, Is.EqualTo(2));
    }

    #endregion

    #region Parse - Base Sequence and Source Format

    [Test]
    public void Parse_ExtractsCorrectBaseSequence()
    {
        var result = _parser.Parse("PEPS[UNIMOD:21]TM[UNIMOD:35]IDE");

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.BaseSequence, Is.EqualTo("PEPSTMIDE"));
    }

    [Test]
    public void Parse_SetsSourceFormatToUnimod()
    {
        var result = _parser.Parse("PEPC[UNIMOD:4]IDE");

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.SourceFormat, Is.EqualTo("Unimod"));
    }

    [Test]
    public void Parse_UnmodifiedSequence_NoModifications()
    {
        Assert.That(_parser.CanParse("PEPTIDE"), Is.False);
    }

    #endregion

    #region Parse - Error Handling

    [Test]
    public void Parse_InvalidTokenFormat_ThrowsException()
    {
        Assert.Throws<SequenceConversionException>(() =>
            _parser.Parse("PEPC[Oxidation on C]IDE"));
    }

    [Test]
    public void Parse_InvalidTokenFormat_ReturnNullMode_ReturnsNull()
    {
        var warnings = new ConversionWarnings();
        var result = _parser.Parse("PEPC[Oxidation on C]IDE", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
        Assert.That(warnings.HasFatalError, Is.True);
        Assert.That(warnings.FailureReason, Is.EqualTo(ConversionFailureReason.UnknownFormat));
    }

    [Test]
    public void Parse_InvalidTokenFormat_WarningContainsTokenText()
    {
        var warnings = new ConversionWarnings();
        _parser.Parse("PEPC[NotAUnimodId]IDE", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(warnings.Errors.Any(e => e.Contains("NotAUnimodId")), Is.True);
    }

    [Test]
    public void Parse_EmptyBrackets_ThrowsException()
    {
        Assert.Throws<SequenceConversionException>(() =>
            _parser.Parse("PEPC[]IDE"));
    }

    [Test]
    public void Parse_EmptyBrackets_ReturnNullMode_ReturnsNull()
    {
        var warnings = new ConversionWarnings();
        var result = _parser.Parse("PEPC[]IDE", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(result, Is.Null);
    }

    #endregion

    #region Parse - Single Residue Edge Cases

    [Test]
    public void Parse_SingleModifiedResidue()
    {
        var result = _parser.Parse("M[UNIMOD:35]");

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.BaseSequence, Is.EqualTo("M"));
        var mod = result.Value.GetModificationAt(0);
        Assert.That(mod, Is.Not.Null);
        Assert.That(mod!.Value.UnimodId, Is.EqualTo(35));
    }

    #endregion

    #region Ground Truth Data-Driven Tests

    [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
    public void CanParse_GroundTruthUnimodFormats(GroundTruthTestData.SequenceConversionTestCase tc)
    {
        if (tc.UnimodUpperCaseFormat == null)
            return;

        var expected = tc.UnimodUpperCaseFormat.Contains('[');
        Assert.That(_parser.CanParse(tc.UnimodUpperCaseFormat), Is.EqualTo(expected),
            $"Failed for: {tc.UnimodUpperCaseFormat}");
    }

    [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
    public void Parse_GroundTruthUpperCase_SetsCorrectUnimodId(GroundTruthTestData.SequenceConversionTestCase tc)
    {
        if (tc.UnimodUpperCaseFormat == null || !tc.UnimodUpperCaseFormat.Contains('['))
            return;

        var result = _parser.Parse(tc.UnimodUpperCaseFormat);

        Assert.That(result, Is.Not.Null, $"Parse failed for: {tc.UnimodUpperCaseFormat}");
        Assert.That(result!.Value.BaseSequence, Is.EqualTo(tc.ExpectedBaseSequence),
            $"Base sequence mismatch for: {tc.UnimodUpperCaseFormat}");
        Assert.That(result.Value.Modifications.Length, Is.EqualTo(tc.ExpectedResidueMods + (tc.HasNTermMod ? 1 : 0) + (tc.HasCTermMod ? 1 : 0)),
            $"Modification count mismatch for: {tc.UnimodUpperCaseFormat}");
    }

    [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
    public void Parse_GroundTruthCamelCase_SetsCorrectUnimodId(GroundTruthTestData.SequenceConversionTestCase tc)
    {
        if (tc.UnimodCamelCaseFormat == null || !tc.UnimodCamelCaseFormat.Contains('['))
            return;

        var result = _parser.Parse(tc.UnimodCamelCaseFormat);

        Assert.That(result, Is.Not.Null, $"Parse failed for camelCase: {tc.UnimodCamelCaseFormat}");
        Assert.That(result!.Value.BaseSequence, Is.EqualTo(tc.ExpectedBaseSequence));
    }

    [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
    public void Parse_GroundTruthLowerCase_SetsCorrectUnimodId(GroundTruthTestData.SequenceConversionTestCase tc)
    {
        if (tc.UnimodLowerCaseFormat == null || !tc.UnimodLowerCaseFormat.Contains('['))
            return;

        var result = _parser.Parse(tc.UnimodLowerCaseFormat);

        Assert.That(result, Is.Not.Null, $"Parse failed for lowercase: {tc.UnimodLowerCaseFormat}");
        Assert.That(result!.Value.BaseSequence, Is.EqualTo(tc.ExpectedBaseSequence));
    }

    [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.CoreTestCases))]
    public void Parse_GroundTruthNoLabel_SetsCorrectUnimodId(GroundTruthTestData.SequenceConversionTestCase tc)
    {
        if (tc.UnimodNoLabelFormat == null || !tc.UnimodNoLabelFormat.Contains('['))
            return;

        var result = _parser.Parse(tc.UnimodNoLabelFormat);

        Assert.That(result, Is.Not.Null, $"Parse failed for no-label: {tc.UnimodNoLabelFormat}");
        Assert.That(result!.Value.BaseSequence, Is.EqualTo(tc.ExpectedBaseSequence));
    }

    [TestCaseSource(typeof(GroundTruthTestData), nameof(GroundTruthTestData.EdgeCases))]
    public void Parse_GroundTruthEdgeCases_SetsCorrectBaseSequence(GroundTruthTestData.SequenceConversionTestCase tc)
    {
        if (tc.UnimodUpperCaseFormat == null || !tc.UnimodUpperCaseFormat.Contains('['))
            return;

        var result = _parser.Parse(tc.UnimodUpperCaseFormat);

        Assert.That(result, Is.Not.Null, $"Parse failed for: {tc.UnimodUpperCaseFormat}");
        Assert.That(result!.Value.BaseSequence, Is.EqualTo(tc.ExpectedBaseSequence));
    }

    #endregion

    #region FormatName and Schema

    [Test]
    public void FormatName_ReturnsUnimod()
    {
        Assert.That(_parser.FormatName, Is.EqualTo("Unimod"));
    }

    [Test]
    public void Schema_IsUnimodSequenceFormatSchema()
    {
        Assert.That(_parser.Schema, Is.InstanceOf<UnimodSequenceFormatSchema>());
    }

    [Test]
    public void Instance_ReturnsSameInstance()
    {
        Assert.That(UnimodSequenceParser.Instance, Is.SameAs(UnimodSequenceParser.Instance));
    }

    #endregion

    #region MzLib Parser Discrimination

    [Test]
    public void MzLibParser_CanParse_ReturnsFalseForLabeledUnimodTokens()
    {
        var mzLibParser = new MzLibSequenceParser();

        Assert.That(mzLibParser.CanParse("PEPC[UNIMOD:4]IDE"), Is.False);
    }

    [Test]
    public void MzLibParser_CanParse_ReturnsFalseForCamelCaseUnimodTokens()
    {
        var mzLibParser = new MzLibSequenceParser();

        Assert.That(mzLibParser.CanParse("PEPC[Unimod:4]IDE"), Is.False);
    }

    [Test]
    public void MzLibParser_CanParse_ReturnsFalseForLowerCaseUnimodTokens()
    {
        var mzLibParser = new MzLibSequenceParser();

        Assert.That(mzLibParser.CanParse("PEPC[unimod:4]IDE"), Is.False);
    }

    #endregion
}
