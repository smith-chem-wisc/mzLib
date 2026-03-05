using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
using Omics;
using Omics.Modifications;
using Omics.Modifications.Conversion;
using Proteomics.ProteolyticDigestion;

namespace Test.Omics;

[TestFixture]
[ExcludeFromCodeCoverage]
public class SequenceConverterTests
{
    [Test]
    public void SequenceConverter_ConvertsMetaMorpheusToUnimod()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        var converted = SequenceConverter.Default.ConvertFullSequence(
            peptide,
            ModificationNamingConvention.Unimod);

        Assert.That(converted, Does.Contain("[Unimod:Carbamidomethyl on C]"));
    }

    [Test]
    public void SequenceConverter_RemovesIncompatibleMods_WhenConfigured()
    {
        ModificationMotif.TryGetMotif("M", out var motif);
        var customMod = new Modification(
            _originalId: "CustomNull",
            _modificationType: "Custom",
            _target: motif,
            _monoisotopicMass: null,
            _chemicalFormula: null);

        Mods.AddOrUpdateModification(customMod);

        var customSequence = $"M[Custom:{customMod.IdWithMotif}]A";
        var peptide = new PeptideWithSetModifications(customSequence, Mods.AllModsKnownDictionary);

        var success = SequenceConverter.Default.TryConvertFullSequence(
            peptide,
            ModificationNamingConvention.Unimod,
            SequenceConversionHandlingMode.RemoveIncompatibleMods,
            out var converted,
            out var reason);

        Assert.That(success, Is.True);
        Assert.That(converted, Is.EqualTo(peptide.BaseSequence));
        Assert.That(reason, Is.EqualTo(SequenceConversionFailureReason.ModificationsRemoved));
    }

    [Test]
    public void SequenceConverter_MassShiftFallback_ReturnsMassShiftSequence()
    {
        ModificationMotif.TryGetMotif("M", out var motif);
        var customMod = new Modification(
            _originalId: "CustomMassShift",
            _modificationType: "CustomMassShift",
            _target: motif,
            _monoisotopicMass: 42.010565,
            _chemicalFormula: null);

        Mods.AddOrUpdateModification(customMod);

        var sequence = $"M[{customMod.ModificationType}:{customMod.IdWithMotif}]A";
        var peptide = new PeptideWithSetModifications(sequence, Mods.AllModsKnownDictionary);

        var success = SequenceConverter.Default.TryConvertFullSequence(
            peptide,
            ModificationNamingConvention.Unimod,
            SequenceConversionHandlingMode.UseMassShifts,
            out var converted,
            out var reason);

        Assert.That(success, Is.True);
        Assert.That(reason, Is.EqualTo(SequenceConversionFailureReason.UsedMassShifts));
        Assert.That(converted, Is.EqualTo(peptide.FullSequenceWithMassShift()));
    }

    [Test]
    public void FullSequenceWithMassShift_DefaultPrecisionUsesFourDecimals()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        Assert.That(peptide.FullSequenceWithMassShifts, Does.Contain("[+57.0215]"));
    }

    [Test]
    public void FullSequenceWithMassShift_AllowsUnsignedFormatting()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        var unsigned = peptide.FullSequenceWithMassShift(decimalPlaces: 2, signed: false);

        Assert.That(unsigned, Does.Contain("[57.02]"));
    }

    [Test]
    public void SequenceConverter_ToMassShiftNotation_UsesConfiguredPrecision()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        var formatted = SequenceConverter.ToMassShiftNotation(peptide.FullSequence);

        Assert.That(formatted, Is.EqualTo(peptide.FullSequenceWithMassShifts));
    }

    [Test]
    public void SequenceConverter_FromMassShiftNotation_ReturnsMatchingModDictionary()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        var massShift = peptide.FullSequenceWithMassShifts;
        var modifications = SequenceConverter.FromMassShiftNotation(massShift);

        Assert.That(modifications.ContainsKey(3), Is.True);
        Assert.That(modifications[3].IdWithMotif, Is.EqualTo("Carbamidomethyl on C"));
    }

    [Test]
    public void SequenceConverter_StringConversion_UsesMassShiftFallback()
    {
        ModificationMotif.TryGetMotif("S", out var motif);
        var customMod = new Modification(
            _originalId: "CustomMassShiftString",
            _modificationType: "CustomMassShiftString",
            _target: motif,
            _monoisotopicMass: 79.966331,
            _chemicalFormula: null);

        Mods.AddOrUpdateModification(customMod);

        var peptide = new PeptideWithSetModifications(
            $"MS[{customMod.ModificationType}:{customMod.IdWithMotif}]T",
            Mods.AllModsKnownDictionary);

        var converted = SequenceConverter.Default.ConvertFullSequence(
            peptide.FullSequence,
            ModificationNamingConvention.Mixed,
            ModificationNamingConvention.UniProt,
            SequenceConversionHandlingMode.UseMassShifts);

        Assert.That(converted, Is.EqualTo(peptide.FullSequenceWithMassShift()));
    }

    [Test]
    public void SequenceConverter_StringConversion_RoundTripsMetaMorpheus()
    {
        var peptide = new PeptideWithSetModifications(
            "PEPTC[Common Fixed:Carbamidomethyl on C]IDE",
            Mods.AllKnownProteinModsDictionary);

        var fullSequence = peptide.FullSequence;

        var same = SequenceConverter.Default.ConvertFullSequence(
            fullSequence,
            ModificationNamingConvention.MetaMorpheus,
            ModificationNamingConvention.MetaMorpheus);

        Assert.That(same, Is.EqualTo(fullSequence));

        var unimod = SequenceConverter.Default.ConvertFullSequence(
            fullSequence,
            ModificationNamingConvention.MetaMorpheus,
            ModificationNamingConvention.Unimod,
            SequenceConversionHandlingMode.UsePrimarySequence);

        Assert.That(unimod, Is.Not.Null);
    }

    [Test]
    public void SequenceConverter_Interface_ConvertsModDictionary()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        ISequenceConverter converter = SequenceConverter.Default;

        var success = converter.TryConvertModifications(
            peptide.AllModsOneIsNterminus,
            peptide.BaseSequence,
            ModificationNamingConvention.MetaMorpheus,
            SequenceConversionHandlingMode.ThrowException,
            out var convertedMods,
            out var reason);

        Assert.That(success, Is.True);
        Assert.That(reason, Is.EqualTo(SequenceConversionFailureReason.None));
        Assert.That(convertedMods, Is.Not.Null);
        Assert.That(convertedMods!.Count, Is.EqualTo(peptide.AllModsOneIsNterminus.Count));
    }

    [Test]
    public void SequenceConverter_PublicConvertModifications_ReturnsDictionary()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        var converted = SequenceConverter.Default.ConvertModifications(
            peptide.AllModsOneIsNterminus,
            peptide.BaseSequence,
            ModificationNamingConvention.Unimod);

        Assert.That(converted.Count, Is.EqualTo(peptide.AllModsOneIsNterminus.Count));
    }

    [Test]
    public void SequenceTargetConverter_ReturnsChronologerSequence()
    {
        var peptide = new PeptideWithSetModifications("PEPTIDE", new Dictionary<string, Modification>());

        var success = SequenceTargetConverter.TryConvert(
            peptide,
            SequenceConversionTarget.Chronologer,
            SequenceConversionHandlingMode.RemoveIncompatibleMods,
            out var converted,
            out var reason);

        Assert.That(success, Is.True);
        Assert.That(converted, Is.EqualTo("-PEPTIDE_"));
        Assert.That(reason, Is.Null);
    }

    [Test]
    public void SequenceTargetConverter_PrimarySequenceFallback()
    {
        ModificationMotif.TryGetMotif("M", out var motif);
        var customMod = new Modification(
            _originalId: "CustomNull",
            _modificationType: "CustomNull",
            _target: motif,
            _monoisotopicMass: null,
            _chemicalFormula: null);

        Mods.AddOrUpdateModification(customMod);

        var peptide = new PeptideWithSetModifications(
            $"M[{customMod.ModificationType}:{customMod.IdWithMotif}]A",
            Mods.AllModsKnownDictionary);

        var success = SequenceTargetConverter.TryConvert(
            peptide,
            SequenceConversionTarget.Chronologer,
            SequenceConversionHandlingMode.UsePrimarySequence,
            out var converted,
            out var reason);

        Assert.That(success, Is.True);
        Assert.That(converted, Is.EqualTo("-MA_"));
        Assert.That(reason, Is.EqualTo(SequenceConversionFailureReason.UsedPrimarySequence));
    }

    [Test]
    public void SequenceConverterExtensions_ReturnMetaMorpheusSequence()
    {
        var peptide = new PeptideWithSetModifications(
            "AC[Common Fixed:Carbamidomethyl on C]DE",
            Mods.AllKnownProteinModsDictionary);

        var converted = peptide.ToMetaMorpheusSequence();

        Assert.That(converted, Is.EqualTo(peptide.FullSequence));
    }
}
