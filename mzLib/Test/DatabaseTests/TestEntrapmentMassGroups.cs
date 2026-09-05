using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases.EntrapmentGeneration;

namespace Test.DatabaseTests;

/// <summary>
/// The invariant in the units a search engine actually sees: per (peptide, MASS GROUP), not per
/// (peptide, modification name). The distinction is the whole point -- "the companion carries the
/// same modifications" and "the companion has the same mass groups" are different sentences, and
/// only the second is what the engine explores.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class EntrapmentMassGroupTests
{
    private static readonly HashSet<string> NothingForbidden = new();

    private static IDigestionParams Tryptic =>
        new DigestionParams("trypsin", minPeptideLength: 7, maxMissedCleavages: 0);

    private static Modification Mod(string motifResidue, string id, double mass,
        string restriction = "Anywhere.")
    {
        ModificationMotif.TryGetMotif(motifResidue, out ModificationMotif motif);
        return new Modification(_originalId: id, _modificationType: "Common Biological",
            _target: motif, _locationRestriction: restriction, _monoisotopicMass: mass);
    }

    // +15.994915 on three different residues. A search sees one mass and three candidate residues;
    // a name-based grouping would see three unrelated modifications.
    private static Modification OxidationOnM() => Mod("M", "Oxidation", 15.994915);
    private static Modification OxidationOnW() => Mod("W", "Oxidation", 15.994915);
    private static Modification HydroxylationOnP() => Mod("P", "Hydroxylation", 15.994915);

    [Test]
    public void OneMassGathersEveryModificationOfThatMassAndAllTheirResidues()
    {
        var index = new MassGroupIndex(new[] { OxidationOnM(), OxidationOnW(), HydroxylationOnP() });

        Assert.That(index.Groups.Count, Is.EqualTo(1),
            "one mass is one group however many names carry it");
        MassGroup group = index.Groups[0];
        Assert.That(group.Residues.OrderBy(c => c), Is.EqualTo(new[] { 'M', 'P', 'W' }));
        Assert.That(group.Modifications.Count, Is.EqualTo(3),
            "three distinct IdWithMotif values -- the names differ, the mass does not, and it is the "
            + "mass a search competes on");
        Assert.That(group.Modifications,
            Is.EquivalentTo(new[] { "Oxidation on M", "Oxidation on W", "Hydroxylation on P" }));
        Assert.That(group.Span, Is.Zero);
    }

    [Test]
    public void TheToleranceDecidesWhetherNearIsobarsAreOneCandidateSetOrTwo()
    {
        // Phospho and sulfo are 9.5 mDa apart -- one candidate set or two depending on what the
        // instrument and the search can separate, which is the consumer's call and not ours.
        var phospho = Mod("S", "Phosphorylation", 79.966331);
        var sulfo = Mod("Y", "Sulfonation", 79.956815);

        Assert.That(new MassGroupIndex(new[] { phospho, sulfo }, 0.0).Groups.Count, Is.EqualTo(2));
        Assert.That(new MassGroupIndex(new[] { phospho, sulfo }, 0.001).Groups.Count, Is.EqualTo(2));
        Assert.That(new MassGroupIndex(new[] { phospho, sulfo }, 0.010).Groups.Count, Is.EqualTo(1),
            "at 10 mDa the two become one candidate set");
    }

    [Test]
    public void AWildcardMotifCanCarryItsMassAnywhere()
    {
        // "Acetylation on X" is kept in real configurations because the dictionary has no
        // residue-specific N-terminal acetyl record. Treating X as a literal residue would make
        // every peptide's capacity for it zero.
        var index = new MassGroupIndex(new[] { Mod("X", "Carbamyl", 43.005814) });

        Assert.That(MassGroupIndex.Capacity("PEPTIDEK", index.Groups[0]), Is.EqualTo(8));
    }

    [Test]
    public void AModificationWithoutAMassDefinesNoGroup()
    {
        ModificationMotif.TryGetMotif("S", out ModificationMotif motif);
        var massless = new Modification(_originalId: "Unknown", _modificationType: "Common Biological",
            _target: motif, _locationRestriction: "Anywhere.");

        // Dropped rather than treated as zero: a zero-mass group would swallow every other
        // unmassed modification and read as a real candidate set.
        Assert.That(new MassGroupIndex(new[] { massless }).Groups, Is.Empty);
    }

    [Test]
    public void ANegativeToleranceIsRefused()
    {
        Assert.Throws<MzLibException>(() => new MassGroupIndex(new[] { OxidationOnM() }, -0.001));
    }

    [Test]
    public void CapacityAndAnnotatedSitesBothSurviveTheRearrangement()
    {
        // The invariant, both clauses, on a protein whose annotations are ordinary "Anywhere." ones.
        const string sequence = "MAAALGGDRSMGVDTTPFAWENDRQITTLGGYK";
        var index = new MassGroupIndex(new[] { OxidationOnM(), OxidationOnW(), HydroxylationOnP() });
        var target = new Protein(sequence, "P12345",
            oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 11, new List<Modification> { OxidationOnM() } },   // the M in the middle piece
                { 20, new List<Modification> { OxidationOnW() } },   // the W in the middle piece
            });
        Assert.That(sequence[10], Is.EqualTo('M'), "fixture: annotation 1 must sit on the M");
        Assert.That(sequence[19], Is.EqualTo('W'), "fixture: annotation 2 must sit on the W");

        Protein companion = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden,
            out EntrapmentAssembly assembly);

        var comparison = new MassGroupComparison(index);
        comparison.Add(target, companion, assembly, Tryptic.MinLength);

        Assert.That(comparison.PeptidesCompared, Is.GreaterThan(0), "the fixture must compare something");
        MassGroupTally tally = comparison.Tallies.Single();
        Assert.That(tally.TargetCapacity, Is.EqualTo(tally.CompanionCapacity));
        Assert.That(tally.TargetAnnotatedSites, Is.EqualTo(tally.CompanionAnnotatedSites));
        Assert.That(tally.TargetAnnotatedSites, Is.EqualTo(2), "both annotations must be counted");
        Assert.That(tally.PeptidesWithCapacityMismatch, Is.Zero);
        Assert.That(tally.PeptidesWithAnnotatedSiteMismatch, Is.Zero);
        Assert.That(comparison.Holds, Is.True);
    }

    [Test]
    public void APeptideTerminalAnnotationIsCountedOnlyWhereASearchCouldApplyIt()
    {
        // The case the whole section exists for, and the one a modification tally cannot see. The
        // annotation is transported and written either way; what decides whether it is a hypothesis
        // is whether it still sits at the peptide's first residue.
        //
        // With the piece-terminal anchor in place this passes; it is written against the SEARCH's
        // view rather than the protein's so that it starts failing again the moment that anchor
        // stops holding, which is the only way this class of defect is ever visible.
        const string sequence = "MAAALGGDRSGGVDTTPFAWENDRQITTLGGYK";
        var waterLoss = Mod("S", "Water Loss", -18.010565, "Peptide N-terminal.");
        var index = new MassGroupIndex(new[] { waterLoss });
        var target = new Protein(sequence, "P12345",
            oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 10, new List<Modification> { waterLoss } }
            });
        Assert.That(sequence[9], Is.EqualTo('S'), "fixture: the annotation opens an interior piece");

        Protein companion = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden,
            out EntrapmentAssembly assembly);

        var comparison = new MassGroupComparison(index);
        comparison.Add(target, companion, assembly, Tryptic.MinLength);

        MassGroupTally tally = comparison.Tallies.Single();
        Assert.That(tally.TargetAnnotatedSites, Is.EqualTo(1),
            "the fixture must produce a usable target site, or it asserts nothing");
        Assert.That(tally.CompanionAnnotatedSites, Is.EqualTo(1),
            "the companion must offer the same hypothesis, not merely carry the same annotation");
        Assert.That(tally.PeptidesWithAnnotatedSiteMismatch, Is.Zero);
        Assert.That(comparison.Holds, Is.True);
    }

    [Test]
    public void AnExcisedPeptideIsNotCountedAsAMatchOrAsAMismatch()
    {
        // A per-peptide equality is satisfied vacuously by a peptide that no longer exists, which is
        // why the totals are published beside it. An excised piece has no companion to compare with,
        // so it contributes to neither -- the achieved fold ratio is where that loss is reported.
        const string sequence = "SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR";
        var index = new MassGroupIndex(new[] { Mod("S", "Phosphorylation", 79.966331) });
        var target = new Protein(sequence, "P12345");

        Protein companion = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden,
            out EntrapmentAssembly assembly);

        Assert.That(assembly.ExcisedCount, Is.EqualTo(1), "fixture must actually excise a piece");
        var comparison = new MassGroupComparison(index);
        comparison.Add(target, companion, assembly, Tryptic.MinLength);

        int searchable = assembly.Pieces.Count(p => p.Outcome != PieceOutcome.Excised
                                                    && p.TargetPiece.Length >= Tryptic.MinLength);
        Assert.That(comparison.PeptidesCompared, Is.EqualTo(searchable));
        Assert.That(comparison.Holds, Is.True);
    }

    [Test]
    public void TheDistributionIsReportedGroupByGroupEvenWhereNothingIsAnnotated()
    {
        // "This mass had no annotated site on either side" and "this mass was not looked at" are
        // different statements, and only a row can tell them apart.
        var index = new MassGroupIndex(new[]
        {
            OxidationOnM(), Mod("S", "Phosphorylation", 79.966331), Mod("K", "Acetylation", 42.010565)
        });
        var target = new Protein("MAAALGGDRSGGVDTTPFAWENDRQITTLGGYK", "P12345");
        Protein companion = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden,
            out EntrapmentAssembly assembly);

        var comparison = new MassGroupComparison(index);
        comparison.Add(target, companion, assembly, Tryptic.MinLength);
        string tsv = comparison.ToTabSeparated();

        Assert.That(comparison.Tallies.Count, Is.EqualTo(3));
        Assert.That(tsv, Does.Contain("# massGroupToleranceDaltons\t0"));
        Assert.That(tsv, Does.Contain("# invariantHolds\ttrue"));
        Assert.That(tsv.Split('\n').Count(l => l.StartsWith("15.99")), Is.EqualTo(1));
        Assert.That(tsv.Split('\n').Count(l => l.StartsWith("42.01")), Is.EqualTo(1));
        Assert.That(tsv.Split('\n').Count(l => l.StartsWith("79.96")), Is.EqualTo(1));
        Assert.That(comparison.Tallies.Select(t => t.Group.RepresentativeMass),
            Is.Ordered, "rows ascend by mass so two reports can be diffed");
    }

    [Test]
    public void TheComparisonRefusesIncompleteInput()
    {
        var index = new MassGroupIndex(new[] { OxidationOnM() });
        var comparison = new MassGroupComparison(index);
        var target = new Protein("MAAALGGDRSGGVDTTPFAWENDRQITTLGGYK", "P12345");
        Protein companion = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden,
            out EntrapmentAssembly assembly);

        Assert.Throws<MzLibException>(() => comparison.Add(null, companion, assembly, 7));
        Assert.Throws<MzLibException>(() => comparison.Add(target, null, assembly, 7));
        Assert.Throws<MzLibException>(() => comparison.Add(target, companion, null, 7));
        Assert.Throws<MzLibException>(() => new MassGroupComparison(null));
    }
}
