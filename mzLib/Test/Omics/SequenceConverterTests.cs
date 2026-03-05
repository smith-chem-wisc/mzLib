using System.Diagnostics.CodeAnalysis;
using NUnit.Framework;
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
            out _);

        Assert.That(success, Is.True);
        Assert.That(converted, Is.EqualTo(peptide.BaseSequence));
    }
}
