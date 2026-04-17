using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
[System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
public class KnownModsLookupTests
{
    #region Constructor Tests

    [Test]
    public void Constructor_NullDictionary_FallsBackToAllKnownMods()
    {
        var lookup = new TestableKnownModsLookup(null!);

        Assert.That(lookup.ExposedCandidateSet, Is.SameAs(Mods.AllKnownMods));
    }

    [Test]
    public void Constructor_EmptyDictionary_ProducesEmptyCandidateSet()
    {
        var lookup = new TestableKnownModsLookup(new Dictionary<string, Modification>());

        Assert.That(lookup.ExposedCandidateSet, Is.Empty);
    }

    [Test]
    public void Constructor_PopulatedDictionary_MaterializesValues()
    {
        var mod = CreateModification("TestMod", 'S');
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };

        var lookup = new TestableKnownModsLookup(dict);

        Assert.That(lookup.ExposedCandidateSet.Count, Is.EqualTo(1));
        Assert.That(lookup.ExposedCandidateSet, Does.Contain(mod));
    }

    [Test]
    public void Constructor_CustomMassTolerance_IsRespected()
    {
        var mod = CreateModification("Tolerant", 'K', monoisotopicMass: 100.0);
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };

        var tight = new KnownModsLookup(dict, massTolerance: 0.001);
        var loose = new KnownModsLookup(dict, massTolerance: 1.0);

        var slightlyOff = CanonicalModification.AtResidue(0, 'K', "[+100.01]", mass: 100.01);

        Assert.That(tight.TryResolve(slightlyOff), Is.Null);
        Assert.That(loose.TryResolve(slightlyOff), Is.Not.Null);
    }

    #endregion

    #region Name Tests

    [Test]
    public void Name_ReturnsKnownMods()
    {
        var lookup = new KnownModsLookup(new Dictionary<string, Modification>());
        Assert.That(lookup.Name, Is.EqualTo("Known Mods"));
    }

    #endregion

    #region Resolution by Name / OriginalId Tests

    [Test]
    public void TryResolve_ByOriginalId_ResolvesFromDictionary()
    {
        var mod = CreateModification("MyMod", 'S');
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'S', "MyMod");
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_ByIdWithMotif_ResolvesFromDictionary()
    {
        var mod = CreateModification("Acetyl", 'K');
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'K', mod.IdWithMotif, mzLibId: mod.IdWithMotif);
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_ByTypePrefixId_ResolvesFromDictionary()
    {
        var mod = CreateModification("Phospho", 'T', modificationType: "Common Variable");
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'T', $"Common Variable:{mod.IdWithMotif}");
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_NameNotInDictionary_ReturnsNull()
    {
        var mod = CreateModification("Known", 'S');
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'S', "Unknown");
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Null);
    }

    #endregion

    #region Resolution by UNIMOD ID Tests

    [Test]
    public void TryResolve_ByUnimodId_ResolvesFromDictionary()
    {
        var motif = GetMotif('C');
        var mod = new Modification(
            _originalId: "Carbamidomethyl",
            _accession: "UNIMOD:4",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("C2H3NO"),
            _databaseReference: new Dictionary<string, IList<string>>
            {
                ["UNIMOD"] = new List<string> { "4" }
            });
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'C', "UNIMOD:4", unimodId: 4);
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_ByUnimodIdInAccession_ResolvesFromDictionary()
    {
        var motif = GetMotif('S');
        var mod = new Modification(
            _originalId: "Oxidation",
            _accession: "UNIMOD:35",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("O"));
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'S', "UNIMOD:35", unimodId: 35);
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_UnimodIdNotInDictionary_ReturnsNull()
    {
        var mod = CreateModification("SomeMod", 'S');
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'S', "UNIMOD:999", unimodId: 999);
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Null);
    }

    #endregion

    #region Resolution by Formula Tests

    [Test]
    public void TryResolve_ByFormula_ResolvesFromDictionary()
    {
        var formula = ChemicalFormula.ParseFormula("O");
        var mod = CreateModification("OxMod", 'M', chemicalFormula: formula);
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'M', "unknown", formula: formula);
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    #endregion

    #region Resolution by Mass Tests

    [Test]
    public void TryResolve_ByMass_ResolvesFromDictionary()
    {
        var mod = CreateModification("Heavy", 'K', monoisotopicMass: 42.0101);
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict, massTolerance: 0.01);

        var canonical = CanonicalModification.AtResidue(0, 'K', "[+42.01]", mass: 42.01);
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_ByMass_OutsideTolerance_ReturnsNull()
    {
        var mod = CreateModification("Precise", 'K', monoisotopicMass: 42.0101);
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict, massTolerance: 0.001);

        var canonical = CanonicalModification.AtResidue(0, 'K', "[+50.0]", mass: 50.0);
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Null);
    }

    #endregion

    #region Terminal Modification Tests

    [Test]
    public void TryResolve_NTerminalModification_ResolvesFromDictionary()
    {
        var mod = new Modification(
            _originalId: "Acetyl",
            _modificationType: "Test",
            _target: GetMotif('X'),
            _locationRestriction: "Peptide N-terminal.");
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtNTerminus("Acetyl on X");
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_CTerminalModification_ResolvesFromDictionary()
    {
        var mod = new Modification(
            _originalId: "Amidation",
            _modificationType: "Test",
            _target: GetMotif('X'),
            _locationRestriction: "Peptide C-terminal.");
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtCTerminus("Amidation on X");
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    #endregion

    #region Already-Resolved Tests

    [Test]
    public void TryResolve_AlreadyResolved_ReturnsUnchanged()
    {
        var mod = CreateModification("Resolved", 'S', modificationType: "Known Mods");
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification
            .AtResidue(0, 'S', "Resolved")
            .WithResolvedModification(mod, 0, ModificationPositionType.Residue);

        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod));
    }

    #endregion

    #region Residue Disambiguation Tests

    [Test]
    public void TryResolve_AmbiguousMultiResidue_ReturnsNull()
    {
        var ser = CreateModification("Shared", 'S');
        var thr = CreateModification("Shared", 'T');
        var dict = new Dictionary<string, Modification>
        {
            { ser.IdWithMotif, ser },
            { thr.IdWithMotif, thr }
        };
        var lookup = new KnownModsLookup(dict);

        var result = lookup.TryResolve("Shared", null);

        Assert.That(result, Is.Null);
    }

    [Test]
    public void TryResolve_ResidueDisambiguatesCorrectly()
    {
        var ser = CreateModification("Shared", 'S');
        var thr = CreateModification("Shared", 'T');
        var dict = new Dictionary<string, Modification>
        {
            { ser.IdWithMotif, ser },
            { thr.IdWithMotif, thr }
        };
        var lookup = new KnownModsLookup(dict);

        var canonicalS = CanonicalModification.AtResidue(0, 'S', "Shared");
        var canonicalT = CanonicalModification.AtResidue(0, 'T', "Shared");

        var resultS = lookup.TryResolve(canonicalS);
        var resultT = lookup.TryResolve(canonicalT);

        Assert.That(resultS, Is.Not.Null);
        Assert.That(resultS!.Value.MzLibModification, Is.SameAs(ser));
        Assert.That(resultT, Is.Not.Null);
        Assert.That(resultT!.Value.MzLibModification, Is.SameAs(thr));
    }

    #endregion

    #region Empty Candidate Set Tests

    [Test]
    public void TryResolve_EmptyCandidateSet_ReturnsNull()
    {
        var lookup = new KnownModsLookup(new Dictionary<string, Modification>());
        var canonical = CanonicalModification.AtResidue(0, 'S', "Anything");

        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Null);
    }

    [Test]
    public void TryResolve_StringOverload_EmptyCandidateSet_ReturnsNull()
    {
        var lookup = new KnownModsLookup(new Dictionary<string, Modification>());

        var result = lookup.TryResolve("Anything", 'S');

        Assert.That(result, Is.Null);
    }

    #endregion

    #region NormalizeRepresentation Tests

    [Test]
    public void NormalizeRepresentation_DoesNotAlterInput()
    {
        var lookup = new KnownModsLookup(new Dictionary<string, Modification>());
        var canonical = CanonicalModification.AtResidue(0, 'S', "  MyMod  ");

        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Null);
    }

    #endregion

    #region Dictionary Independence Tests

    [Test]
    public void Constructor_DoesNotMutateSourceDictionary()
    {
        var mod = CreateModification("Safe", 'S');
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var originalCount = dict.Count;

        var lookup = new KnownModsLookup(dict);

        Assert.That(dict.Count, Is.EqualTo(originalCount));
    }

    [Test]
    public void Constructor_ReflectsLiveViewNotSnapshot()
    {
        var mod = CreateModification("Snapshot", 'S');
        var dict = new Dictionary<string, Modification> { { mod.IdWithMotif, mod } };
        var lookup = new KnownModsLookup(dict);

        var mod2 = CreateModification("AddedLater", 'T');
        dict.Add(mod2.IdWithMotif, mod2);

        var canonical = CanonicalModification.AtResidue(0, 'T', "AddedLater");
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null, "Dictionary.ValueCollection is a live view; additions are visible.");
        Assert.That(result!.Value.MzLibModification, Is.SameAs(mod2));
    }

    #endregion

    #region Multiple Mods Same Name Different Type Tests

    [Test]
    public void TryResolve_TypePrefixDisambiguatesSameNameDifferentResidue()
    {
        var modS = CreateModification("MyMod", 'S', modificationType: "TypeA");
        var modT = CreateModification("MyMod", 'T', modificationType: "TypeA");
        var dict = new Dictionary<string, Modification>
        {
            { modS.IdWithMotif, modS },
            { modT.IdWithMotif, modT }
        };
        var lookup = new KnownModsLookup(dict);

        var canonical = CanonicalModification.AtResidue(0, 'S', $"TypeA:{modS.IdWithMotif}", mzLibId: $"TypeA:{modS.IdWithMotif}");
        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.Not.Null);
        Assert.That(result!.Value.MzLibModification, Is.SameAs(modS));
    }

    #endregion

    private class TestableKnownModsLookup : KnownModsLookup
    {
        public TestableKnownModsLookup(Dictionary<string, Modification> knownMods, double massTolerance = 0.001)
            : base(knownMods, massTolerance)
        {
        }

        public IReadOnlyCollection<Modification> ExposedCandidateSet => CandidateSet;
    }

    private static Modification CreateModification(
        string originalId,
        char residue,
        string location = "Anywhere.",
        double? monoisotopicMass = 42.0101,
        string modificationType = "Test",
        ChemicalFormula? chemicalFormula = null)
    {
        return new Modification(
            _originalId: originalId,
            _modificationType: modificationType,
            _target: GetMotif(residue),
            _locationRestriction: location,
            _chemicalFormula: chemicalFormula,
            _monoisotopicMass: monoisotopicMass);
    }

    private static ModificationMotif GetMotif(char residue)
    {
        if (!ModificationMotif.TryGetMotif(residue.ToString(), out var motif))
        {
            throw new InvalidOperationException($"Unable to create motif for '{residue}'.");
        }

        return motif;
    }
}
