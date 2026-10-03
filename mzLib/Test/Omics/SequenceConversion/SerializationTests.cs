using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.ProForma;
using System;
using System.Collections.Generic;
using System.Linq;
using Transcriptomics;
using Transcriptomics.Digestion;
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

    #region MzLib Serializer Output Read Back By PeptideWithSetModifications

    private static ISequenceParser GetParser(string parser) =>
        parser == "ProForma" ? ProFormaSequenceParser.Instance : MzLibSequenceParser.Instance;

    private static Modification KnownMod(string modificationType, string idWithMotif, string locationRestriction) =>
        Mods.AllKnownMods.First(m => m.ModificationType == modificationType && m.IdWithMotif == idWithMotif && m.LocationRestriction == locationRestriction);

    /// <summary>
    /// The peptide with <paramref name="baseSequence"/> that digesting <paramref name="protein"/> produces with
    /// <paramref name="variableMod"/> (if any) on it.
    /// </summary>
    private static PeptideWithSetModifications Digested(string protein, string baseSequence, Modification? variableMod, string protease = "trypsin") =>
        new Protein(protein, "P")
            .Digest(new DigestionParams(protease: protease, maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                new List<Modification>(), variableMod == null ? new List<Modification>() : new List<Modification> { variableMod })
            .Single(p => p.BaseSequence == baseSequence && p.AllModsOneIsNterminus.Count > 0);

    /// <summary>
    /// The peptide as one of mzLib's writers puts it, read back by the parser for that writer's format.
    /// </summary>
    private static CanonicalSequence Written(PeptideWithSetModifications peptide, string writer) => writer switch
    {
        "ProFormaWriter" => ProFormaSequenceParser.Instance.Parse(ProFormaSequenceSerializer.Instance.Serialize(peptide.ToCanonicalSequence())!)!.Value,
        "MassShifts" => MassShiftSequenceParser.Instance.Parse(peptide.FullSequenceWithMassShifts)!.Value,
        "FullSequence" => MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value,
        _ => peptide.ToCanonicalSequence()
    };

    [Test]
    [TestCase("[UNIMOD:737]-PEPTIDEK", "ProForma", 1, 229.1629)]
    [TestCase("PEPM[UNIMOD:35]IDEK", "ProForma", 5, 15.9949)]
    public void MzLibSerializer_ModificationMzLibCannotReadBack_IsResolvedToOneItCan(string input, string parser, int oneIsNTerminusIndex, double expectedMass)
    {
        var canonical = GetParser(parser).Parse(input);
        Assert.That(canonical, Is.Not.Null);

        var result = MzLibSequenceSerializer.Instance.Serialize(canonical.Value);

        Assert.That(result, Is.Not.Null);
        var peptide = new PeptideWithSetModifications(result);
        Assert.That(peptide.BaseSequence, Is.EqualTo(canonical.Value.BaseSequence));
        Assert.That(peptide.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { oneIsNTerminusIndex }));
        Assert.That(peptide.AllModsOneIsNterminus[oneIsNTerminusIndex].MonoisotopicMass, Is.EqualTo(expectedMass).Within(0.001));
    }

    [Test]
    [TestCase("[Unimod:TMT6plex on X]PEPTIDEK")]
    [TestCase("PEPM[Common Variable:Oxidation on M]IDEK")]
    public void MzLibSerializer_NameMzLibCanReadBack_IsWrittenUnchanged(string input)
    {
        var canonical = MzLibSequenceParser.Instance.Parse(input);
        Assert.That(canonical, Is.Not.Null);

        var result = MzLibSequenceSerializer.Instance.Serialize(canonical.Value);

        Assert.That(result, Is.EqualTo(input));
        Assert.That(() => new PeptideWithSetModifications(result), Throws.Nothing);
    }

    [Test]
    public void MzLibSerializer_MzLibNameNoLookupKnows_FailsPerHandlingMode()
    {
        var canonical = MzLibSequenceParser.Instance.Parse("PEPN[N-Glycosylation:H5N2 on N]K");
        Assert.That(canonical, Is.Not.Null);

        AssertFailsPerHandlingMode(canonical.Value, "PEPNK");
    }

    // Each modification is written by a real writer. What the lookup resolves it to from that text, or the
    // dictionary entry its written name reads as, has no name that reads back as it, or would be read back at
    // another position.
    [Test]
    [TestCase("ProFormaWriter", "Less Common", "Methylation on X", "N-terminal.", "AGGGGK", "AGGGGK")]
    [TestCase("ProFormaWriter", "Less Common", "Ethylation on X", "Peptide N-terminal.", "AGGGGK", "AGGGGK")]
    [TestCase("MassShifts", "Less Common", "Ethylation on X", "Peptide N-terminal.", "AGGGGK", "AGGGGK")]
    [TestCase("MassShifts", "UniProt", "N,N-dimethylalanine on A", "N-terminal.", "AGGGGK", "AGGGGK")]
    [TestCase("MassShifts", "UniProt", "N2,N2-dimethylarginine on R", "N-terminal.", "RGGGGK", "R")]
    [TestCase("MassShifts", "UniProt", "Lysine methyl ester on K", "C-terminal.", "GGGK", "GGGK")]
    [TestCase("ToCanonicalSequence", "Less Common", "Methylation on X", "C-terminal.", "GGGA", "GGGA")]
    [TestCase("FullSequence", "Less Common", "Methylation on X", "C-terminal.", "GGGA", "GGGA")]
    public void MzLibSerializer_WrittenModificationWithNoNameThatReadsBackAsIt_FailsPerHandlingMode(string writer, string modificationType, string idWithMotif, string locationRestriction, string protein, string baseSequence)
    {
        var peptide = Digested(protein, baseSequence, KnownMod(modificationType, idWithMotif, locationRestriction));

        AssertFailsPerHandlingMode(Written(peptide, writer), baseSequence);
    }

    [Test]
    public void MzLibSerializer_CnbrNTerminalLactoneWrittenAsMassShift_FailsPerHandlingMode()
    {
        // The lactone's mass on the first residue resolves to a C-terminal entry, which reading would move to the end.
        var peptide = Digested("AAAMPEPTIDE", "MPEPTIDE", null, "CNBr_N");

        AssertFailsPerHandlingMode(Written(peptide, "MassShifts"), "MPEPTIDE");
    }

    private static void AssertFailsPerHandlingMode(CanonicalSequence canonical, string expectedWithoutIt, params int[] remainingOneIsNTerminusIndices)
    {
        Assert.That(() => MzLibSequenceSerializer.Instance.Serialize(canonical, null, SequenceConversionHandlingMode.ThrowException),
            Throws.TypeOf<SequenceConversionException>());
        Assert.That(MzLibSequenceSerializer.Instance.Serialize(canonical, null, SequenceConversionHandlingMode.ReturnNull), Is.Null);

        var warnings = new ConversionWarnings();
        var removed = MzLibSequenceSerializer.Instance.Serialize(canonical, warnings, SequenceConversionHandlingMode.RemoveIncompatibleElements);
        Assert.That(removed, Is.EqualTo(expectedWithoutIt));
        Assert.That(warnings.HasIncompatibleItems, Is.True);
        Assert.That(new PeptideWithSetModifications(removed).AllModsOneIsNterminus.Keys, Is.EquivalentTo(remainingOneIsNTerminusIndices));
    }

    [Test]
    [TestCase(SequenceConversionHandlingMode.ThrowException)]
    [TestCase(SequenceConversionHandlingMode.ReturnNull)]
    [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements)]
    public void MzLibSerializer_CnbrDigestedPeptide_ReadsBackWithItsHomoserineLactoneAtTheCTerminus(SequenceConversionHandlingMode mode)
    {
        // Digestion puts the protease's C-terminal homoserine lactone on the last residue, not at the C-terminus.
        var homoserineLactone = ProteaseDictionary.Dictionary["CNBr"].CleavageMod;
        var peptide = Digested("AAAMPEPTIDEMKKK", "PEPTIDEM", null, "CNBr");
        Assert.That(peptide.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 9 }));

        var result = MzLibSequenceSerializer.Instance.Serialize(peptide.ToCanonicalSequence(), null, mode);

        Assert.That(result, Is.EqualTo("PEPTIDEM[Protease:Homoserine lactone on M]"));
        var readBack = new PeptideWithSetModifications(result,
            new Dictionary<string, Modification> { { homoserineLactone.IdWithMotif, homoserineLactone } });
        Assert.That(readBack.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 10 }));
        Assert.That(readBack.AllModsOneIsNterminus[10].MonoisotopicMass, Is.EqualTo(homoserineLactone.MonoisotopicMass).Within(1e-5));
    }

    [Test]
    public void MzLibSerializer_CnbrPeptideWrittenAsProForma_ConvertsBackWithItsHomoserineLactoneAtTheCTerminus()
    {
        // The ProForma writer puts the lactone on the last residue as its UNIMOD id; the lookup resolves that to a
        // C-terminal entry, which must stay writable there.
        var homoserineLactone = ProteaseDictionary.Dictionary["CNBr"].CleavageMod;
        var peptide = Digested("AAAMPEPTIDEMKKK", "PEPTIDEM", null, "CNBr");
        var proForma = ProFormaSequenceSerializer.Instance.Serialize(peptide.ToCanonicalSequence());
        Assert.That(proForma, Is.EqualTo("PEPTIDEM[UNIMOD:11]"));

        var result = MzLibSequenceSerializer.Instance.Serialize(ProFormaSequenceParser.Instance.Parse(proForma)!.Value, null, SequenceConversionHandlingMode.ThrowException);

        var readBack = new PeptideWithSetModifications(result);
        Assert.That(readBack.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 10 }));
        Assert.That(readBack.AllModsOneIsNterminus[10].MonoisotopicMass, Is.EqualTo(homoserineLactone.MonoisotopicMass).Within(1e-5));
    }

    // The default instance resolves an oligo's masses among RNA modifications too, so mzLib's oligo writer output
    // reads back as the same oligo.
    [Test]
    [TestCase("Biological", "2'-O-Methyladenosine on A", "GUAACUG")]
    [TestCase("Metal", "Sodium on A", "GUAACUG")]
    public void MzLibSerializer_OligoWrittenWithMassShifts_ReadsBackAsTheSameOligo(string modificationType, string idWithMotif, string sequence)
    {
        var mod = Mods.MetaMorpheusRnaModifications.First(m => m.ModificationType == modificationType && m.IdWithMotif == idWithMotif);
        var oligo = new RNA(sequence)
            .Digest(new RnaDigestionParams(maxMods: 1), new List<Modification>(), new List<Modification> { mod })
            .First(o => o.AllModsOneIsNterminus.Count > 0);
        var digested = oligo.AllModsOneIsNterminus.Single();

        var result = MzLibSequenceSerializer.Instance.Serialize(MassShiftSequenceParser.Instance.Parse(global::Omics.BioPolymerWithSetModsExtensions.FullSequenceWithMassShift(oligo))!.Value,
            null, SequenceConversionHandlingMode.ThrowException);

        var readBack = new OligoWithSetMods(result);
        Assert.That(readBack.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { digested.Key }));
        Assert.That(readBack.AllModsOneIsNterminus[digested.Key].MonoisotopicMass, Is.EqualTo(digested.Value.MonoisotopicMass).Within(1e-3));
    }

    [Test]
    [TestCase(SequenceConversionHandlingMode.ThrowException)]
    [TestCase(SequenceConversionHandlingMode.ReturnNull)]
    [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements)]
    public void MzLibSerializer_AttachedModificationMzLibDoesNotKnow_IsWrittenByNameForAReaderThatHasIt(SequenceConversionHandlingMode mode)
    {
        ModificationMotif.TryGetMotif("N", out var motif);
        var glycan = new Modification(_originalId: "H5N2", _modificationType: "N-Glycosylation", _target: motif,
            _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("C46H76N4O35"));
        Assert.That(Mods.AllKnownProteinModsDictionary.ContainsKey(glycan.IdWithMotif) || Mods.AllKnownRnaModsDictionary.ContainsKey(glycan.IdWithMotif), Is.False);
        var dictionary = new Dictionary<string, Modification> { { glycan.IdWithMotif, glycan } };
        const string expected = "PEPN[N-Glycosylation:H5N2 on N]K";
        var peptide = new PeptideWithSetModifications(expected, dictionary);

        var result = MzLibSequenceSerializer.Instance.Serialize(peptide.ToCanonicalSequence(), null, mode);

        Assert.That(result, Is.EqualTo(expected));
        Assert.That(new PeptideWithSetModifications(result, dictionary).AllModsOneIsNterminus[5], Is.SameAs(glycan));
        Assert.That(() => new PeptideWithSetModifications(result), Throws.TypeOf<MzLibUtil.MzLibException>());
    }

    [Test]
    public void MzLibSerializer_AttachedModificationEquivalentToTheDictionaryEntry_IsWrittenByName()
    {
        // MetaMorpheus's own oxidation is not the dictionary's (UNIMOD's) entry for its name, but reads back with
        // the same mass and terminus. A lookup with no candidates leaves the attached modification as the only answer.
        var metaMorpheusOxidation = KnownMod("Common Variable", "Oxidation on M", "Anywhere.");
        Assert.That(metaMorpheusOxidation, Is.Not.SameAs(Mods.AllKnownProteinModsDictionary["Oxidation on M"]));
        var peptide = Digested("PEPMIDEK", "PEPMIDEK", metaMorpheusOxidation);
        var serializer = new MzLibSequenceSerializer(new UnimodModificationLookup(Array.Empty<Modification>()));

        var result = serializer.Serialize(peptide.ToCanonicalSequence(), null, SequenceConversionHandlingMode.ThrowException);

        Assert.That(result, Is.EqualTo("PEPM[Common Variable:Oxidation on M]IDEK"));
    }

    // The ProForma writers give these the UNIMOD id from their database reference, whose UNIMOD entry has another
    // mass, so the ProForma text itself names a different modification.
    private static readonly HashSet<string> ProFormaAccessionWithAnotherMass = new()
    {
        "N6,N6,N6-trimethyl-5-hydroxylysine on K",
        "N-linked (Lac) (glycation) lysine on K",
    };

    [Test]
    public void MzLibSerializer_EveryProteinModification_WrittenAsProForma_ConvertsToAnMzLibSequenceThatReadsBackAsIt()
    {
        // Each protein modification is put on a peptide by digestion and written by the ProForma writer. A conversion
        // either fails per the handling mode or reads back where and as heavy as the digested modification (a
        // C-terminal one digestion put on the last residue reads back at the C-terminus).
        var converted = 0;
        Assert.Multiple(() =>
        {
            foreach (var mod in Mods.AllProteinModsList.Where(m => m.Target != null && m.MonoisotopicMass.HasValue && m.ValidModification
                                                                   && !ProFormaAccessionWithAnotherMass.Contains(m.IdWithMotif)))
            {
                var residues = string.Concat(mod.Target.ToString().Select(c => c == 'X' || c == 'x' ? 'A' : char.ToUpperInvariant(c)));
                var peptides = new Protein(residues + "GGGGK" + "GGG" + residues, "P")
                    .Digest(new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                        new List<Modification>(), new List<Modification> { mod })
                    .Where(p => p.AllModsOneIsNterminus.Count > 0);
                foreach (var peptide in peptides)
                {
                    var proForma = ProFormaSequenceSerializer.Instance.Serialize(peptide.ToCanonicalSequence(), null, SequenceConversionHandlingMode.ReturnNull);
                    if (proForma == null)
                        continue;

                    // The writer can name a modification in a way its own parser rejects; that is the writer's to fix.
                    var parsed = ProFormaSequenceParser.Instance.Parse(proForma, null, SequenceConversionHandlingMode.ReturnNull);
                    if (parsed == null)
                        continue;

                    var result = MzLibSequenceSerializer.Instance.Serialize(parsed.Value, null, SequenceConversionHandlingMode.ReturnNull);
                    if (result == null)
                        continue;

                    var length = peptide.BaseSequence.Length;
                    var expected = peptide.AllModsOneIsNterminus.ToDictionary(
                        kvp => kvp.Key == length + 1 && kvp.Value.LocationRestriction.Contains("C-terminal") ? length + 2 : kvp.Key,
                        kvp => kvp.Value.MonoisotopicMass!.Value);
                    PeptideWithSetModifications? readBack = null;
                    Assert.That(() => readBack = new PeptideWithSetModifications(result), Throws.Nothing, $"{peptide.FullSequence} -> {proForma} -> {result}");
                    if (readBack == null)
                        continue;

                    Assert.That(readBack.AllModsOneIsNterminus.Keys, Is.EquivalentTo(expected.Keys), $"{peptide.FullSequence} -> {proForma} -> {result}");
                    foreach (var (index, mass) in expected)
                    {
                        if (readBack.AllModsOneIsNterminus.TryGetValue(index, out var readBackMod))
                            Assert.That(readBackMod.MonoisotopicMass, Is.EqualTo(mass).Within(0.01), $"{peptide.FullSequence} -> {proForma} -> {result}");
                    }
                    converted++;
                }
            }
        });
        Assert.That(converted, Is.GreaterThan(1000));
    }

    #endregion
}
