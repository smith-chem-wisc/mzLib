using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using MassSpectrometry;
using Omics.Modifications;
using Omics.SequenceConversion;
using Omics.Fragmentation;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class ModomicsSequenceConversionTests
{
    private static readonly ModomicsSequenceParser Parser = ModomicsSequenceParser.Instance;

    [TestCase(null, false)]
    [TestCase("", false)]
    [TestCase("   ", false)]
    [TestCase("GUACUG", false)]
    [TestCase("GJACUG", true)]
    [TestCase(" GJ ACUG ", true)]
    [TestCase("GU\u2603AC", false)]
    [TestCase("P", false)]
    public void CanParse_RecognizesModomicsAlphabetAndCodes(string input, bool expected)
    {
        Assert.That(Parser.CanParse(input), Is.EqualTo(expected));
    }

    [Test]
    public void CanParse_RejectsUnknownCharactersAfterValidSequence()
    {
        Assert.That(Parser.CanParse("GJACUG\u2603"), Is.False);
    }

    [TestCase("GUACUG", "GUACUG", 0)]
    [TestCase("GJACUGCBUCUA#UGAA#CA", "GUACUGCCUCUAGUGAAGCA", 4)]
    [TestCase("/UCCAGU#CAGUACJG", "AUCCAGUGCAGUACUG", 3)]
    [TestCase("UUCAAGUA:UCCAGGAUAGGCU", "UUCAAGUAAUCCAGGAUAGGCU", 1)]
    [TestCase("UUCAAGUA=UCCAGGAUAGGCU", "UUCAAGUAAUCCAGGAUAGGCU", 1)]
    [TestCase("GA:C", "GAAC", 1)]
    [TestCase("GA=C", "GAAC", 1)]
    [TestCase("G:C", "GAC", 1)]
    [TestCase("[G]", "AGU", 2)]
    public void ParsesWorkbookSequences(string input, string expectedBaseSequence, int expectedModificationCount)
    {
        var sequence = Parser.Parse(input);

        Assert.That(sequence, Is.Not.Null);
        Assert.That(sequence!.Value.BaseSequence, Is.EqualTo(expectedBaseSequence));
        Assert.That(sequence.Value.ModificationCount, Is.EqualTo(expectedModificationCount));
        Assert.That(sequence.Value.AllModificationsResolved, Is.True);
    }

    [Test]
    public void FullSequenceParser_TrimsWhitespaceAfterModomicsNamespace()
    {
        var oligo = new OligoWithSetMods("A[Modomics: 2'-O-methyladenosine on A]U");

        Assert.That(oligo.BaseSequence, Is.EqualTo("AU"));
        Assert.That(oligo.AllModsOneIsNterminus, Does.ContainKey(2));
        Assert.That(oligo.AllModsOneIsNterminus[2].OriginalId,
            Is.EqualTo("2'-O-methyladenosine"));
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
    public void Parse_WhitespaceOnlyInput_ReturnsNullInReturnNullMode()
    {
        var warnings = new ConversionWarnings();

        var sequence = Parser.Parse(" \t\r\n", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(sequence, Is.Null);
        Assert.That(warnings.HasFatalError, Is.True);
    }

    [Test]
    public void Parse_WhitespaceOnlyInput_ThrowsInThrowExceptionMode()
    {
        Assert.Throws<SequenceConversionException>(() => Parser.Parse(" \t\r\n"));
    }

    [Test]
    public void Parse_FivePrimeTerminalCodeAtSequenceStart_CreatesNTerminalModification()
    {
        var terminalCode = Mods.ModomicsLoadReport.ModificationsByAbbreviation
            .First(pair => pair.Value.Any(modification =>
                Mods.ModomicsLoadReport.TerminalModifications.Contains(modification)))
            .Key[0];

        var sequence = Parser.Parse($"{terminalCode}AC");

        Assert.That(sequence, Is.Not.Null);
        Assert.That(sequence!.Value.BaseSequence, Is.EqualTo("AC"));
        Assert.That(sequence.Value.Modifications.Length, Is.EqualTo(1));
        Assert.That(sequence.Value.Modifications[0].PositionType, Is.EqualTo(ModificationPositionType.NTerminus));
        Assert.That(sequence.Value.Modifications[0].OriginalRepresentation, Is.EqualTo(terminalCode.ToString()));
    }

    [Test]
    public void Parse_FivePrimeTerminalCodeAfterResidue_ReturnsNullInReturnNullMode()
    {
        var terminalCode = Mods.ModomicsLoadReport.ModificationsByAbbreviation
            .First(pair => pair.Value.Any(modification =>
                Mods.ModomicsLoadReport.TerminalModifications.Contains(modification)))
            .Key[0];
        var warnings = new ConversionWarnings();

        var sequence = Parser.Parse($"A{terminalCode}C", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(sequence, Is.Null);
        Assert.That(warnings.HasFatalError, Is.True);
        Assert.That(warnings.IncompatibleItems, Does.Contain(terminalCode.ToString()));
    }

    [Test]
    public void Parse_AmbiguousCodeWithUnmatchedTarget_ReturnsNullInReturnNullMode()
    {
        var ambiguousCode = Mods.ModomicsLoadReport.ModificationsByAbbreviation
            .FirstOrDefault(pair => pair.Key.Length == 1
                && !"ACGUPY".Contains(pair.Key[0])
                && pair.Value.Count > 1)
            .Key;

        Assume.That(ambiguousCode, Is.Not.Null.And.Not.Empty);
        var warnings = new ConversionWarnings();

        var sequence = Parser.Parse($"A{ambiguousCode}Z", warnings, SequenceConversionHandlingMode.ReturnNull);

        Assert.That(sequence, Is.Null);
        Assert.That(warnings.HasFatalError, Is.True);
        Assert.That(warnings.IncompatibleItems, Does.Contain(ambiguousCode));
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
    public void BracketCodesResolveToTheirOwnTargetResidues()
    {
        var leftBracket = ModomicsModificationLookup.Instance.TryResolve("[", 'A');
        var rightBracket = ModomicsModificationLookup.Instance.TryResolve("]", 'U');

        Assert.That(leftBracket, Is.Not.Null);
        Assert.That(leftBracket!.Value.MzLibModification!.OriginalId,
            Is.EqualTo("2-methylthio-N6-threonylcarbamoyladenosine"));
        Assert.That(rightBracket, Is.Not.Null);
        Assert.That(rightBracket!.Value.MzLibModification!.OriginalId,
            Is.EqualTo("1-methylpseudouridine"));
    }

    [Test]
    public void EveryLoadedOneLetterCodeResolvesThroughModomicsLookup()
    {
        foreach (var abbreviation in Mods.ModomicsLoadReport.ModificationsByAbbreviation.Keys)
        {
            Assert.That(abbreviation.Length, Is.EqualTo(1), $"Unexpected multi-character code: {abbreviation}");
            foreach (var candidate in Mods.ModomicsLoadReport.ModificationsByAbbreviation[abbreviation])
            {
                var target = candidate.Target.Motif[0];
                Assert.That(ModomicsModificationLookup.Instance.TryResolveCode(abbreviation[0], target, out var modification),
                    Is.True, $"Code {abbreviation} did not resolve for target {target}");
                Assert.That(modification!.ModificationType, Is.AnyOf("Modomics", "5' Terminal Cap"));
            }
        }
    }

    [Test]
    public void ServiceRegistersModomicsSourceAndConverter()
    {
        Assert.That(SequenceConversionService.Default.AvailableSourceFormats,
            Does.Contain("Modomics"));
        Assert.That(SequenceConversionService.Default.AvailableConverters,
            Does.Contain("Modomics-mzLib"));
    }

    [TestCase(
        "GJACUGCBUCUA#UGAA#CA",
         "GU[Modomics:2'-O-methyluridine on U]ACUGCC[Modomics:2'-O-methylcytidine on C]UCUAG[Modomics:2'-O-methylguanosine on G]UGAAG[Modomics:2'-O-methylguanosine on G]CA")]
    [TestCase(
        "/UCCAGU#CAGUACJG",
         "A[Modomics:2-methyladenosine on A]UCCAGUG[Modomics:2'-O-methylguanosine on G]CAGUACU[Modomics:2'-O-methyluridine on U]G")]
    [TestCase(
        "UUCAAGUA:UCCAGGAUAGGCU",
        "UUCAAGUAA[Biological: 2'-O-Methyladenosine on A]UCCAGGAUAGGCU")]
     [TestCase(
         "UUCAAGUA=UCCAGGAUAGGCU",
         "UUCAAGUAA[Modomics:N6-methyladenosine on A]UCCAGGAUAGGCU")]
    [TestCase(
        "UCCCUGAGACCCUA:CUUGUGA",
        "UCCCUGAGACCCUAA[Common Biological: Methylation on A]CUUGUGA")]
    public void ModomicsAndMetaMorpheusConstructionProduceIdenticalFragments(
        string modomicsSequence,
        string metaMorpheusFullSequence)
    {
        var modomicsCanonical = SequenceConversionService.Default.Parse(modomicsSequence, "Modomics");
        var modomicsToMzlib = SequenceConversionService.Default.Convert(modomicsSequence, "Modomics", "mzLib");

        Assert.That(modomicsCanonical, Is.Not.Null);

        var metaMorpheusOligo = new OligoWithSetMods(metaMorpheusFullSequence);
        var convertedOligo = new OligoWithSetMods(modomicsToMzlib);

        Assert.That(convertedOligo.BaseSequence, Is.EqualTo(metaMorpheusOligo.BaseSequence));
        Assert.That(convertedOligo.AllModsOneIsNterminus.Keys, Is.EquivalentTo(metaMorpheusOligo.AllModsOneIsNterminus.Keys));

        foreach (var dissociationType in new[] { DissociationType.CID, DissociationType.HCD })
        {
            var metaMorpheusProducts = new List<Product>();
            var convertedProducts = new List<Product>();
            metaMorpheusOligo.Fragment(dissociationType, FragmentationTerminus.Both, metaMorpheusProducts);
            convertedOligo.Fragment(dissociationType, FragmentationTerminus.Both, convertedProducts);

            Assert.That(convertedProducts.Count, Is.EqualTo(metaMorpheusProducts.Count), dissociationType.ToString());
            var expected = metaMorpheusProducts.OrderBy(product => product.ProductType)
                .ThenBy(product => product.Terminus)
                .ThenBy(product => product.FragmentNumber)
                .ThenBy(product => product.NeutralLoss)
                .ToList();
            var actual = convertedProducts.OrderBy(product => product.ProductType)
                .ThenBy(product => product.Terminus)
                .ThenBy(product => product.FragmentNumber)
                .ThenBy(product => product.NeutralLoss)
                .ToList();

            for (var i = 0; i < expected.Count; i++)
            {
                Assert.That(actual[i].ProductType, Is.EqualTo(expected[i].ProductType), dissociationType.ToString());
                Assert.That(actual[i].Terminus, Is.EqualTo(expected[i].Terminus), dissociationType.ToString());
                Assert.That(actual[i].FragmentNumber, Is.EqualTo(expected[i].FragmentNumber), dissociationType.ToString());
                Assert.That(actual[i].NeutralLoss, Is.EqualTo(expected[i].NeutralLoss).Within(1e-9), dissociationType.ToString());
                Assert.That(actual[i].NeutralMass, Is.EqualTo(expected[i].NeutralMass).Within(1e-9), dissociationType.ToString());
            }
        }
    }
}
