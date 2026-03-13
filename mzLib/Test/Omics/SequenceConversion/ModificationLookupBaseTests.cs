using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion;

[TestFixture]
public class ModificationLookupBaseTests
{
    [Test]
    public void TryResolve_ReturnsCanonicalWhenAlreadyResolved()
    {
        var resolved = CreateModification("Resolved", 'S', modificationType: TestLookup.LookupName);
        var canonical = CanonicalModification
            .AtResidue(0, 'S', "Resolved")
            .WithResolvedModification(resolved, 0, ModificationPositionType.Residue);
        var lookup = new TestLookup(new[] { resolved });

        var result = lookup.TryResolve(canonical);

        Assert.That(result, Is.EqualTo(canonical));
    }

    [Test]
    public void TryResolve_UsesMassRepresentationWhenAvailable()
    {
        var oxidation = CreateModification("Oxidation", 'M', monoisotopicMass: 15.9949);
        var lookup = new TestLookup(new[] { oxidation });
        var canonical = CanonicalModification.AtResidue(2, 'M', "[+15.9949]", mass: 15.9949);

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(oxidation));
    }

    [Test]
    public void TryResolve_ReturnsNullWhenNameDoesNotMatchAndNoOtherInfo()
    {
        var known = CreateModification("Known", 'S');
        var lookup = new TestLookup(new[] { known });
        var canonical = CanonicalModification.AtResidue(0, 'S', "Unknown");

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Null);
    }

    [Test]
    public void FilterByName_ExpandsIsobaricAndNTerminalVariants()
    {
        var tmtPro = CreateModification("TMTpro", 'K');
        var tmt18 = CreateModification("TMT18", 'K');
        var lookup = new TestLookup(new[] { tmtPro, tmt18 });

        var matches = lookup.FilterByNamePublic(null, "TMTpro on N-terminus", 'K').ToList();

        Assert.That(matches, Does.Contain(tmtPro));
        Assert.That(matches, Does.Contain(tmt18));
    }

    [Test]
    public void FilterByTerm_RespectsLocationRestrictions()
    {
        var nTerm = CreateModification("NTerm", 'S', location: "N-terminal.");
        var cTerm = CreateModification("CTerm", 'S', location: "C-terminal.");
        var anywhere = CreateModification("Anywhere", 'S', location: "Anywhere.");
        var lookup = new TestLookup(new[] { nTerm, cTerm, anywhere });

        var nTermMatches = lookup.FilterByTermPublic(null, ModificationPositionType.NTerminus).ToList();
        Assert.That(nTermMatches, Does.Contain(nTerm));
        Assert.That(nTermMatches, Does.Contain(anywhere));
        Assert.That(nTermMatches, Does.Not.Contain(cTerm));

        var cTermMatches = lookup.FilterByTermPublic(null, ModificationPositionType.CTerminus).ToList();
        Assert.That(cTermMatches, Does.Contain(cTerm));
        Assert.That(cTermMatches, Does.Contain(anywhere));
        Assert.That(cTermMatches, Does.Not.Contain(nTerm));
    }

    [Test]
    public void FilterByUnimodId_FindsMatchesByAccessionOrDatabaseReference()
    {
        var motif = GetMotif('S');
        var byAccession = new Modification(
            _originalId: "Oxidation",
            _accession: "UNIMOD:35",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("O"));

        var byReference = new Modification(
            _originalId: "Custom",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("C"),
            _databaseReference: new Dictionary<string, IList<string>>
            {
                ["UNIMOD"] = new List<string> { "UNIMOD:36" }
            });

        var lookup = new TestLookup(new[] { byAccession, byReference });

        var accessionMatches = lookup.FilterByUnimodIdPublic(null, 35).ToList();
        Assert.That(accessionMatches, Does.Contain(byAccession));

        var referenceMatches = lookup.FilterByUnimodIdPublic(null, 36).ToList();
        Assert.That(referenceMatches, Does.Contain(byReference));
    }

    [Test]
    public void BuildCacheKey_CombinesAllKeyAttributes()
    {
        var lookup = new TestLookup(Array.Empty<Modification>());
        var formula = ChemicalFormula.ParseFormula("C2H2");
        var key = lookup.BuildCacheKeyPublic(
            "Name",
            'S',
            formula,
            42.010123,
            new AbsoluteTolerance(0.123456),
            ModificationPositionType.NTerminus);

        Assert.That(key, Does.Contain("Name"));
        Assert.That(key, Does.Contain("|res:S"));
        Assert.That(key, Does.Contain("|cf:C2H2"));
        Assert.That(key, Does.Match(".*\\|mass:42\\.0101.*"));
        Assert.That(key, Does.Contain("@tol:0.123456"));
        Assert.That(key, Does.Contain("|term:NTerminus"));
    }

    [Test]
    public void MatchesIdentifier_SupportsTypePrefixesAndSuffixes()
    {
        var mod = CreateModification("Acetyl", 'K', modificationType: "TypeA");
        var lookup = new TestLookup(new[] { mod });

        Assert.That(lookup.MatchesIdentifierPublic(mod, $"TypeA:{mod.IdWithMotif}"), Is.True);
        Assert.That(lookup.MatchesIdentifierPublic(mod, $"prefix:{mod.OriginalId}"), Is.True);
        Assert.That(lookup.MatchesIdentifierPublic(mod, "other:value"), Is.False);
    }

    [Test]
    public void ExpandNameCandidates_ProducesResidueAndIsobaricVariants()
    {
        var lookup = new TestLookup(Array.Empty<Modification>());

        var variants = lookup.ExpandNameCandidatesPublic("TMT18 on N-terminus", 'K').ToList();

        Assert.That(variants, Does.Contain("TMT18 on K"));
        Assert.That(variants, Does.Contain("TMTpro on N-terminus"));
        Assert.That(variants.Any(v => v.Contains("TMTpro on K")), Is.True);
    }

    [Test]
    public void TryResolve_StringPathResolvesToCanonical()
    {
        var formula = ChemicalFormula.ParseFormula("H2O1");
        var mod = CreateModification("WaterLoss", 'S', chemicalFormula: formula);
        var lookup = new TestLookup(new[] { mod });

        var resolved = lookup.TryResolve("WaterLoss on S", 'S', formula, ModificationPositionType.Residue);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibId, Is.EqualTo(mod.IdWithMotif));
    }

    [Test]
    public void TryResolve_FormulaFilterResolvesWhenNameFails()
    {
        var formula = ChemicalFormula.ParseFormula("C2H3");
        var mod = CreateModification("FormulaMatch", 'T', chemicalFormula: formula);
        var lookup = new TestLookup(new[] { mod });
        var canonical = CanonicalModification.AtResidue(2, 'T', "not-a-match", formula: formula);

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_MassFilterUsesTolerance()
    {
        var mod = CreateModification("MassMatch", 'K', monoisotopicMass: 42.02);
        var lookup = new TestLookup(new[] { mod }, massTolerance: 0.05);
        var canonical = CanonicalModification.AtResidue(0, 'K', "[+42.00]", mass: 41.98);

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(mod));
    }

    [Test]
    public void TryResolve_UsesFormulaAndResidueFallbackWhenNameUnknown()
    {
        var formula = ChemicalFormula.ParseFormula("C2H3");
        var ser = CreateModification("Ambiguous", 'S', chemicalFormula: formula);
        var thr = CreateModification("Ambiguous", 'T', chemicalFormula: formula);
        var lookup = new TestLookup(new[] { ser, thr });

        var canonical = CanonicalModification.AtResidue(0, 'T', "UnknownName", formula: formula);

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(thr));
    }

    [Test]
    public void TryResolve_MassFilterCachesNormalizedTolerance()
    {
        var modification = CreateModification("MassCache", 'S', monoisotopicMass: 123.456789);
        var lookup = new CountingMassLookup(new[] { modification }, massTolerance: 0.005);
        var canonical = CanonicalModification.AtResidue(0, 'S', "[+123.456789]", mass: 123.456789);
        var slightlyDifferent = canonical.WithMass(123.456788);

        var firstResolution = lookup.TryResolve(canonical);
        var secondResolution = lookup.TryResolve(slightlyDifferent);

        Assert.That(firstResolution, Is.Not.Null);
        Assert.That(secondResolution, Is.Not.Null);
        Assert.That(firstResolution!.Value.MzLibModification, Is.SameAs(modification));
        Assert.That(secondResolution!.Value.MzLibModification, Is.SameAs(modification));
        Assert.That(lookup.MassFilterCalls, Is.EqualTo(1));
    }

    [Test]
    public void TryResolve_TargetResidueFiltersCandidates()
    {
        var ser = CreateModification("ResidueSpecific", 'S');
        var thr = CreateModification("ResidueSpecific", 'T');
        var lookup = new TestLookup(new[] { ser, thr });
        var canonical = CanonicalModification.AtResidue(0, 'S', "ResidueSpecific");

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(ser));
    }

    [Test]
    public void SelectBestCandidate_FiltersIsobaricWithoutDiagnostics()
    {
        var diag = new Dictionary<DissociationType, List<double>>
        {
            { DissociationType.HCD, new List<double> { 100.1 } }
        };

        var withDiag = new Modification(
            _originalId: "TMTpro",
            _modificationType: "TestDiag",
            _target: GetMotif('K'),
            _locationRestriction: "Anywhere.",
            _diagnosticIons: diag,
            _monoisotopicMass: 304.2);

        var withoutDiag = new Modification(
            _originalId: "TMTpro",
            _modificationType: "Test",
            _target: GetMotif('K'),
            _locationRestriction: "Anywhere.",
            _monoisotopicMass: 304.2);

        var lookup = new TestLookup(new[] { withoutDiag, withDiag });
        var canonical = CanonicalModification.AtResidue(0, 'K', "TMTpro");

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(withDiag));
    }

    [Test]
    public void SelectBestCandidate_PrefersShortestNameWhenChemicallyEquivalent()
    {
        var motif = GetMotif('M');
        var longName = new Modification(
            _originalId: "Methionine (R)-sulfoxide",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("O"));
        var shortName = new Modification(
            _originalId: "Methionine sulfoxide",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("O"));

        var lookup = new TestLookup(new[] { longName, shortName });
        var canonical = CanonicalModification.AtResidue(5, 'M', "Methionine sulfoxide");

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(shortName));
    }

    [Test]
    public void SelectBestCandidate_ReturnsNullWhenAmbiguousResidues()
    {
        var ser = CreateModification("SharedName", 'S');
        var thr = CreateModification("SharedName", 'T');
        var lookup = new TestLookup(new[] { ser, thr });

        var resolved = lookup.TryResolve("SharedName", null);

        Assert.That(resolved, Is.Null);
    }

    [Test]
    public void MatchesTermRestriction_HandlesTerminalStrings()
    {
        var lookup = new TestLookup(Array.Empty<Modification>());

        Assert.That(lookup.MatchesTermRestrictionPublic("Peptide N-terminal.", ModificationPositionType.NTerminus), Is.True);
        Assert.That(lookup.MatchesTermRestrictionPublic("Peptide N-terminal.", ModificationPositionType.CTerminus), Is.False);
        Assert.That(lookup.MatchesTermRestrictionPublic("Peptide C-terminal.", ModificationPositionType.CTerminus), Is.True);
        Assert.That(lookup.MatchesTermRestrictionPublic("Peptide C-terminal.", ModificationPositionType.NTerminus), Is.False);
        Assert.That(lookup.MatchesTermRestrictionPublic("Anywhere.", ModificationPositionType.Residue), Is.True);
    }

    [Test]
    public void IsMassRepresentationDetectsBracketedMasses()
    {
        var lookup = new TestLookup(Array.Empty<Modification>());

        Assert.That(lookup.IsMassRepresentationPublic("[+15.99]"), Is.True);
        Assert.That(lookup.IsMassRepresentationPublic("15.99"), Is.True);
        Assert.That(lookup.IsMassRepresentationPublic("Acetyl"), Is.False);
    }

    [Test]
    public void TryResolve_CachesNullResults()
    {
        var lookup = new CountingTestLookup(new[] { CreateModification("Known", 'S') });

        var first = lookup.TryResolve("Missing", 'S');
        var second = lookup.TryResolve("Missing", 'S');

        Assert.That(first, Is.Null);
        Assert.That(second, Is.Null);
        Assert.That(lookup.FilterByNameCalls, Is.EqualTo(1));
    }

    [Test]
    public void Constructor_MaterializesEnumerables()
    {
        var lookup = new TestLookup(YieldCandidates(CreateModification("Materialized", 'S')));

        Assert.That(lookup.RawCandidateSet.Count, Is.EqualTo(1));
    }

    [Test]
    public void Constructor_DefaultsToKnownModsWhenCandidatesNull()
    {
        var lookup = new TestLookup();

        Assert.That(lookup.RawCandidateSet, Is.SameAs(Mods.AllKnownMods));
    }

    [Test]
    public void FilterByIdentifier_MatchesTypePrefixedNotation()
    {
        var mod = CreateModification("Identified", 'S', modificationType: "CustomType");
        var lookup = new TestLookup(new[] { mod });

        var matches = lookup.FilterByIdentifierPublic(null, $"CustomType:{mod.IdWithMotif}").ToList();

        Assert.That(matches, Does.Contain(mod));
    }

    [Test]
    public void FilterByMotif_IncludesWildcardAndEmptyTargets()
    {
        var specific = CreateModification("Specific", 'S');
        var wildcard = new Modification(
            _originalId: "Wildcard",
            _modificationType: "Test",
            _target: GetMotif('X'),
            _locationRestriction: "Anywhere.");
        var empty = new Modification(
            _originalId: "Empty",
            _modificationType: "Test",
            _target: null,
            _locationRestriction: "Anywhere.");
        var lookup = new TestLookup(new[] { specific, wildcard, empty });

        var matches = lookup.FilterByMotifPublic(new[] { specific, wildcard, empty }, 'S').ToList();

        Assert.That(matches, Does.Contain(specific));
        Assert.That(matches, Does.Contain(wildcard));
        Assert.That(matches, Does.Not.Contain(empty));
    }

    [Test]
    public void GetOverlapScore_ReturnsLongestOverlapLength()
    {
        var lookup = new TestLookup();

        var score = lookup.GetOverlapScorePublic("Acetyl on K", "AcetylK");

        Assert.That(score, Is.EqualTo(6));
    }

    [Test]
    public void GetPrimaryCandidates_UsesUnimodIdWhenAvailable()
    {
        var motif = GetMotif('S');
        var candidate = new Modification(
            _originalId: "Unimod",
            _modificationType: "Test",
            _target: motif,
            _locationRestriction: "Anywhere.",
            _databaseReference: new Dictionary<string, IList<string>>
            {
                ["UNIMOD"] = new List<string> { "UNIMOD:35" }
            });
        var canonical = CanonicalModification.AtResidue(0, 'S', "placeholder", unimodId: 35);
        var lookup = new TestLookup(new[] { candidate });

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(candidate));
    }

    [Test]
    public void TryResolve_ContinuesWhenResolvedToDifferentDatabase()
    {
        var motif = GetMotif('S');
        var lookupMod = CreateModification("Resolved", 'S', modificationType: TestLookup.LookupName);
        var otherDb = new Modification(
            _originalId: "Resolved",
            _modificationType: "OtherDB",
            _target: motif,
            _locationRestriction: "Anywhere.");
        var canonical = CanonicalModification
            .AtResidue(0, 'S', "Resolved")
            .WithResolvedModification(otherDb, 0);
        var lookup = new TestLookup(new[] { lookupMod });

        var resolved = lookup.TryResolve(canonical);

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.Value.MzLibModification, Is.SameAs(lookupMod));
    }

    [Test]
    public void BuildCacheKey_OmitsOptionalPartsWhenMissing()
    {
        var lookup = new TestLookup();
        var key = lookup.BuildCacheKeyPublic("Simple", null, null, null, new AbsoluteTolerance(0.5), null);

        Assert.That(key, Is.EqualTo("Simple"));
    }

    [Test]
    public void FilterByMass_ReturnsEmptyWhenOutsideTolerance()
    {
        var mod = CreateModification("MassSpecific", 'S', monoisotopicMass: 15.99);
        var lookup = new TestLookup(new[] { mod });

        var results = lookup.FilterByMassPublic(new[] { mod }, 100).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void FilterByName_UsesSuffixAfterColon()
    {
        var mod = CreateModification("ColonName", 'S');
        var lookup = new TestLookup(new[] { mod });

        var matches = lookup.FilterByNamePublic(null, $"Type:{mod.IdWithMotif}", null).ToList();

        Assert.That(matches, Does.Contain(mod));
    }

    [Test]
    public void ExpandIsobaricTags_HandlesHyphenatedAndPlainPlex()
    {
        var lookup = new TestLookup();

        var hyphenated = lookup.ExpandNameCandidatesPublic("TMT-plex", null).ToList();
        var plain = lookup.ExpandNameCandidatesPublic("iTRAQplex", null).ToList();

        Assert.That(hyphenated.Any(v => v.Contains("TMTplex", StringComparison.OrdinalIgnoreCase)), Is.True);
        Assert.That(hyphenated.Any(v => v.Equals("TMT", StringComparison.OrdinalIgnoreCase)), Is.True);
        Assert.That(plain.Any(v => v.Contains("-plex", StringComparison.OrdinalIgnoreCase)), Is.True);
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

    private static IEnumerable<Modification> YieldCandidates(params Modification[] modifications)
    {
        foreach (var modification in modifications)
        {
            yield return modification;
        }
    }

    private sealed class CountingTestLookup : TestLookup
    {
        public CountingTestLookup(IEnumerable<Modification> candidates)
            : base(candidates)
        {
        }

        public int FilterByNameCalls { get; private set; }

        protected override IEnumerable<Modification> FilterByName(IEnumerable<Modification> source, string name, char? targetResidue)
        {
            FilterByNameCalls++;
            return base.FilterByName(source, name, targetResidue);
        }
    }

    private sealed class CountingMassLookup : TestLookup
    {
        public CountingMassLookup(IEnumerable<Modification> candidates, double massTolerance)
            : base(candidates, massTolerance: massTolerance)
        {
        }

        public int MassFilterCalls { get; private set; }

        protected override IEnumerable<Modification> FilterByMass(IEnumerable<Modification> source, double mass)
        {
            MassFilterCalls++;
            return base.FilterByMass(source, mass);
        }
    }

    private class TestLookup : ModificationLookupBase
    {
        public const string LookupName = "TestLookup";
        private readonly Func<CanonicalModification, IEnumerable<Modification>>? _primaryFactory;

        public TestLookup(
            IEnumerable<Modification>? candidates = null,
            Func<CanonicalModification, IEnumerable<Modification>>? primaryFactory = null,
            double massTolerance = 0.01)
            : base(candidates, massTolerance)
        {
            _primaryFactory = primaryFactory;
        }

        public override string Name => LookupName;

        public IReadOnlyCollection<Modification> RawCandidateSet => CandidateSet;

        protected override IEnumerable<Modification> GetPrimaryCandidates(CanonicalModification mod)
        {
            if (_primaryFactory != null)
            {
                return _primaryFactory(mod);
            }

            return base.GetPrimaryCandidates(mod);
        }

        public IEnumerable<Modification> FilterByNamePublic(IEnumerable<Modification>? source, string name, char? residue)
            => FilterByName(source ?? CandidateSet, name, residue);

        public IEnumerable<Modification> FilterByMotifPublic(IEnumerable<Modification>? source, char motif)
            => FilterByMotif(source ?? CandidateSet, motif);

        public IEnumerable<Modification> FilterByMassPublic(IEnumerable<Modification>? source, double mass)
            => FilterByMass(source ?? CandidateSet, mass);

        public IEnumerable<Modification> FilterByTermPublic(IEnumerable<Modification>? source, ModificationPositionType term)
            => FilterByTerm(source ?? CandidateSet, term);

        public IEnumerable<Modification> FilterByUnimodIdPublic(IEnumerable<Modification>? source, int unimodId)
            => FilterByUnimodId(source ?? CandidateSet, unimodId);

        public string BuildCacheKeyPublic(string representation, char? residue, ChemicalFormula? formula, double? mass, Tolerance tolerance, ModificationPositionType? term)
            => BuildCacheKey(representation, residue, formula, mass, tolerance, term);

        public bool MatchesIdentifierPublic(Modification modification, string identifier)
            => MatchesIdentifier(modification, identifier);

        public IEnumerable<string> ExpandNameCandidatesPublic(string representation, char? residue)
            => ExpandNameCandidates(representation, residue);

        public IEnumerable<Modification> FilterByIdentifierPublic(IEnumerable<Modification>? source, string identifier)
            => FilterByIdentifier(source ?? CandidateSet, identifier);

        public bool MatchesTermRestrictionPublic(string? restriction, ModificationPositionType term)
            => MatchesTermRestriction(restriction, term);

        public bool IsMassRepresentationPublic(string representation)
            => IsMassRepresentation(representation);

        public int GetOverlapScorePublic(string idWithMotif, string trimmedName)
            => GetOverlapScore(idWithMotif, trimmedName);
    }
}
