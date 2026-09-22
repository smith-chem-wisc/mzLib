using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class ModomicsSequenceConversionTests
{
    private static readonly ModomicsSequenceParser Parser = ModomicsSequenceParser.Instance;

    [TestCase("GUACUG", "GUACUG", 0)]
    [TestCase("GJACUGCBUCUA#UGAA#CA", "GUACUGCCUCUAGUGAAGCA", 4)]
    [TestCase("/UCCAGU#CAGUACJG", "AUCCAGUGCAGUACUG", 3)]
    [TestCase("UUCAAGUA:UCCAGGAUAGGCU", "UUCAAGUAAUCCAGGAUAGGCU", 1)]
    [TestCase("UUCAAGUA=UCCAGGAUAGGCU", "UUCAAGUAAUCCAGGAUAGGCU", 1)]
    public void ParsesWorkbookSequences(string input, string expectedBaseSequence, int expectedModificationCount)
    {
        var sequence = Parser.Parse(input);

        Assert.That(sequence, Is.Not.Null);
        Assert.That(sequence!.Value.BaseSequence, Is.EqualTo(expectedBaseSequence));
        Assert.That(sequence.Value.ModificationCount, Is.EqualTo(expectedModificationCount));
        Assert.That(sequence.Value.AllModificationsResolved, Is.True);
    }

    [Test]
    public void PseudouridineBecomesMzLibYResidue()
    {
        var sequence = Parser.Parse("GPACU");

        Assert.That(sequence, Is.Not.Null);
        Assert.That(sequence!.Value.BaseSequence, Is.EqualTo("GYACU"));
        Assert.That(sequence.Value.ModificationCount, Is.Zero);
    }

    [Test]
    public void ConvertsToTypedMzLibNotation()
    {
        var warnings = new ConversionWarnings();
        var converted = SequenceConversionService.Default.Convert(
            "GJACUGCBUCUA#UGAA#CA",
            "Modomics",
            "mzLib",
            warnings);

        Assert.That(converted, Does.StartWith("GU[Modomics:"));
        Assert.That(converted, Does.Contain("Modomics:"));
        Assert.That(converted, Does.EndWith("]CA"));
        Assert.That(warnings.IsClean, Is.True);
    }

    [Test]
    public void AutoDetectionIdentifiesModomicsButNotPlainRna()
    {
        Assert.That(SequenceConversionService.Default.DetectFormat("GUACUG"), Is.EqualTo("mzLib"));
        Assert.That(SequenceConversionService.Default.DetectFormat("GJACUG"), Is.EqualTo("Modomics"));
    }

    [Test]
    public void UnknownCodeReturnsNullWithoutThrowingInReturnNullMode()
    {
        var warnings = new ConversionWarnings();

        var sequence = Parser.Parse("GU\u2603AC", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(sequence, Is.Null);
        Assert.That(warnings.HasFatalError, Is.True);
        Assert.That(warnings.IncompatibleItems, Does.Contain("\u2603"));
    }

    [Test]
    public void UnknownCodeThrowsInThrowExceptionMode()
    {
        Assert.Throws<SequenceConversionException>(() => Parser.Parse("GU\u2603AC"));
    }

    [Test]
    public void LookupResolvesModomicsInstanceForMethyladenosineCode()
    {
        var resolved = ModomicsModificationLookup.Instance.TryResolve("=", 'A');

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.Not.Null);
        Assert.That(resolved.Value.MzLibModification!.ModificationType, Is.EqualTo("Modomics"));
        Assert.That(resolved.Value.MzLibModification.Target.Motif, Is.EqualTo("A"));
    }

    [Test]
    public void ModomicsCodesAreCaseSensitive()
    {
        Assert.That(ModomicsModificationLookup.Instance.TryResolve("j", 'U'), Is.Null);
        Assert.That(ModomicsModificationLookup.Instance.TryResolve("J", 'U'), Is.Not.Null);
    }

    [Test]
    public void ServiceRegistersModomicsSourceAndConverter()
    {
        Assert.That(SequenceConversionService.Default.AvailableSourceFormats,
            Does.Contain("Modomics"));
        Assert.That(SequenceConversionService.Default.AvailableConverters,
            Does.Contain("Modomics-mzLib"));
    }
}
