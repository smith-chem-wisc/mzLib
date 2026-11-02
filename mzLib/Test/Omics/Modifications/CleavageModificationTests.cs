using System;
using System.Collections.Generic;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;

namespace Test.Omics.Modifications;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class CleavageModificationTests
{
    private static Modification CreateSourceModification()
    {
        ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
        return new Modification(
            _originalId: "Test Cleavage Modification",
            _modificationType: "Digestion Termini",
            _target: motif,
            _locationRestriction: "Oligo 3'-terminal.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H-2 O-1"));
    }

    [TestCase(true, true)]
    [TestCase(false, false)]
    public void Constructor_InvalidFixedVariableCombination_Throws(bool isFixed, bool isVariable)
    {
        Assert.Throws<ArgumentException>(() =>
            new CleavageModification(isFixed, isVariable, CreateSourceModification()));
    }

    [TestCase(true, false)]
    [TestCase(false, true)]
    public void Constructor_ValidFixedVariableCombination_PreservesFlags(bool isFixed, bool isVariable)
    {
        var modification = new CleavageModification(isFixed, isVariable, CreateSourceModification());

        Assert.That(modification.IsFixedMod, Is.EqualTo(isFixed));
        Assert.That(modification.IsVariableMod, Is.EqualTo(isVariable));
        Assert.That(modification.OriginalId, Is.EqualTo("Test Cleavage Modification"));
        Assert.That(modification.LocationRestriction, Is.EqualTo("Oligo 3'-terminal."));
    }

    [Test]
    public void Equality_PlainAndCleavageModifications_AreDistinctTypes()
    {
        var source = CreateSourceModification();
        var cleavageModification = new CleavageModification(true, false, source);

        Assert.That(cleavageModification, Is.Not.EqualTo(source));
        Assert.That(cleavageModification.CompareTo(source), Is.Not.EqualTo(0));

        var modifications = new HashSet<Modification> { source, cleavageModification };
        Assert.That(modifications, Has.Count.EqualTo(2));
    }
}
