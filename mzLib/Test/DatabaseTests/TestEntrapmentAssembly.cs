using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases.EntrapmentGeneration;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class EntrapmentAssemblyTests
{
    private static readonly HashSet<string> NothingForbidden = new();

    private static IDigestionParams Tryptic(int minLength = 7) =>
        new DigestionParams("trypsin", minPeptideLength: minLength, maxMissedCleavages: 2);

    /// <summary>Base pieces of a sequence under the same partition the assembler uses.</summary>
    private static List<string> BasePieces(string sequence, IDigestionParams digestionParams)
    {
        List<int> sites = digestionParams.DigestionAgent.GetDigestionSiteIndices(sequence);
        sites.Sort();
        return Enumerable.Range(0, sites.Count - 1)
            .Select(i => sequence.Substring(sites[i], sites[i + 1] - sites[i]))
            .ToList();
    }

    private static string Sorted(string s) => string.Concat(s.OrderBy(c => c));

    [Test]
    public void Assemble_ProducesOnePermutedPieceForEachTargetPiece()
    {
        const string target = "SYKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK";
        IDigestionParams digestion = Tryptic();

        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, digestion, NothingForbidden);

        List<string> targetPieces = BasePieces(target, digestion);
        List<string> entrapmentPieces = BasePieces(assembly.EntrapmentSequence, digestion);

        Assert.That(assembly.ExcisedCount, Is.Zero);
        Assert.That(entrapmentPieces.Count, Is.EqualTo(targetPieces.Count),
            "the entrapment protein must digest into the same number of pieces");
        for (int i = 0; i < targetPieces.Count; i++)
        {
            Assert.That(Sorted(entrapmentPieces[i]), Is.EqualTo(Sorted(targetPieces[i])),
                "piece " + i + " must be isomeric with its target counterpart");
        }
    }

    [Test]
    public void Assemble_KeepsTheWholeProteinIsomeric()
    {
        const string target = "SYKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK";
        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden);

        Assert.That(assembly.EntrapmentSequence, Is.Not.EqualTo(target));
        Assert.That(Sorted(assembly.EntrapmentSequence), Is.EqualTo(Sorted(target)));
        Assert.That(assembly.EntrapmentSequence.Count(c => c is 'S' or 'T'),
            Is.EqualTo(target.Count(c => c is 'S' or 'T')));
    }

    [Test]
    public void Assemble_KeepsASubLengthPieceVerbatimRatherThanExcisingIt()
    {
        // "AK" has no permutation once K is pinned, but it is below the minimum peptide length, so
        // it is never identifiable on its own. Keeping it verbatim puts no real target peptide into
        // the entrapment search space, and identity is a permutation, so isomerism still holds.
        // Excising every such piece instead would gut the database.
        const string target = "AKSYKALADQMNLLLSKGGVDTTPFAWENDR";
        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden);

        Assert.That(assembly.ExcisedCount, Is.Zero, "a sub-length piece must not be excised");
        Assert.That(assembly.KeptVerbatimCount, Is.EqualTo(1), "exactly the \"AK\" piece");
        Assert.That(assembly.EntrapmentSequence, Does.StartWith("AK"));
        Assert.That(assembly.EntrapmentSequence.Length, Is.EqualTo(target.Length));
    }

    [Test]
    public void Assemble_ExcisesASearchablePieceThatHasNoPermutation()
    {
        // "SSSSSSR" is long enough to be identified and has exactly one arrangement, so there is no
        // honest partner for it. It is removed rather than emitted unchanged, because emitting it
        // would place a genuine target peptide inside the entrapment database.
        const string target = "SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR";
        IDigestionParams digestion = Tryptic();

        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, digestion, NothingForbidden);

        Assert.That(assembly.ExcisedCount, Is.EqualTo(1));
        Assert.That(assembly.EntrapmentSequence, Does.Not.Contain("SSSSSSR"));
        Assert.That(assembly.EntrapmentSequence.Length, Is.EqualTo(target.Length - 7));

        // Every other piece digests exactly as it would have without the excision.
        List<string> targetPieces = BasePieces(target, digestion);
        List<string> entrapmentPieces = BasePieces(assembly.EntrapmentSequence, digestion);
        Assert.That(entrapmentPieces.Count, Is.EqualTo(targetPieces.Count - 1));
    }

    [Test]
    public void Assemble_NeverLeavesATargetPeptideInTheEntrapmentSequence()
    {
        const string target = "SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR";
        IDigestionParams digestion = Tryptic();
        var targetPeptides = BasePieces(target, digestion).ToHashSet();

        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, digestion, targetPeptides);

        foreach (string piece in BasePieces(assembly.EntrapmentSequence, digestion))
        {
            if (piece.Length >= digestion.MinLength)
            {
                Assert.That(targetPeptides, Does.Not.Contain(piece),
                    "'" + piece + "' is a real target peptide and must not appear as entrapment");
            }
        }
    }

    [Test]
    public void Assemble_CountsMissedCleavagePeptidesThatSpanAnExcision()
    {
        // Base-peptide compositions add, so a missed-cleavage peptide is normally isomeric with its
        // target counterpart for free. One that spans an excision is not -- it has no counterpart at
        // all. Those are counted rather than passed off as isomeric, because someone will compute
        // mass statistics over the entrapment set assuming the invariant holds everywhere.
        const string target = "SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR";
        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden);

        // Pieces are SYK | ALADQMNLLLSK | SSSSSSR | GGVDTTPFAWENDR, and the third is excised, so the
        // retained ordinals are 0, 1, 3. At two missed cleavages the runs that straddle that gap are
        // {0,1,3} and {1,3} -- exactly two. An exact count, not "more than none": the arithmetic here
        // is nested and off-by-one prone, and a loose assertion would not notice it drifting.
        Assert.That(assembly.ExcisedCount, Is.EqualTo(1));
        Assert.That(assembly.MissedCleavagePeptidesSpanningAnExcision, Is.EqualTo(2));
    }

    [Test]
    public void Assemble_MapsEveryRetainedPositionAndMarksExcisedOnes()
    {
        const string target = "SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR";
        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden);

        Assert.That(assembly.TargetToEntrapmentPosition.Length, Is.EqualTo(target.Length));

        for (int i = 0; i < target.Length; i++)
        {
            int mapped = assembly.TargetToEntrapmentPosition[i];
            if (mapped < 0)
            {
                continue;   // excised
            }

            Assert.That(assembly.EntrapmentSequence[mapped], Is.EqualTo(target[i]),
                "position " + i + " must map to a position holding the same residue");
        }

        // The excised stretch, and only it, is unmapped.
        Assert.That(assembly.TargetToEntrapmentPosition.Count(p => p < 0), Is.EqualTo(7));
    }

    [Test]
    public void Assemble_RejectsAnEmptySequence()
    {
        var thrown = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentAssembler.Assemble("", Tryptic(), NothingForbidden));
        Assert.That(thrown!.Message, Does.Contain("empty"));
    }

    [Test]
    public void Assemble_CountsBrokenRunsFromTheGapItself()
    {
        // Retained ordinals here are 0, 2, 3 -- the excised piece sits second rather than last.
        // The run {2,3} is contiguous and must NOT count, which is what distinguishes a real
        // difference of ordinals from any other arithmetic that happens to agree on [0,1,3].
        const string target = "ALADQMNLLLSKSSSSSSRGGVDTTPFAWENDRQISTLGGYK";
        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden);

        Assert.That(assembly.ExcisedCount, Is.EqualTo(1));
        Assert.That(assembly.MissedCleavagePeptidesSpanningAnExcision, Is.EqualTo(2));
    }

    [Test]
    public void Assemble_IsDeterministic()
    {
        const string target = "SYKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK";

        string first = EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden).EntrapmentSequence;
        EntrapmentAssembler.Assemble("ACDEFGHIKLMNPQR", Tryptic(), NothingForbidden);  // unrelated work

        Assert.That(EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden).EntrapmentSequence,
            Is.EqualTo(first));
        Assert.That(EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden, seed: 2).EntrapmentSequence,
            Is.Not.EqualTo(first));
    }

    [Test]
    public void Assemble_ReproducesTheCandidateSiteDistributionOfTheTarget()
    {
        // The check per-peptide equality structurally cannot make. Equality is pointwise over
        // peptides present in both databases, so a peptide that was dropped satisfies it vacuously;
        // only the distribution notices an assembly bug. This is the assertion that caught the
        // rule which excised every short piece.
        string[] proteins =
        {
            "SYKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK",
            "AKMTTSTTPAPTTTKQEEEKKKGGVDTTPFAWENDRSSSSSSR",
            "MPSSTSRLLGVDKAAATPTGGGSQSCCCKYPERDNRSTTTK",
            "QLQGASWELQSLRPSLDQLAAHPWMLGADGGVPESCDLRTTK"
        };
        IDigestionParams digestion = Tryptic();

        var expectedSites = new List<int>();
        var entrapmentSites = new List<int>();
        int excised = 0;

        foreach (string protein in proteins)
        {
            EntrapmentAssembly assembly = EntrapmentAssembler.Assemble(protein, digestion, NothingForbidden);

            // Excision is the ONLY thing allowed to change the distribution, so the target's
            // contribution is its searchable pieces minus the ones that had no partner.
            expectedSites.AddRange(assembly.Pieces
                .Where(p => p.Outcome != PieceOutcome.Excised)
                .Select(p => p.TargetPiece)
                .Where(p => p.Length >= digestion.MinLength)
                .Select(p => p.Count(c => c is 'S' or 'T')));
            excised += assembly.ExcisedCount;

            entrapmentSites.AddRange(BasePieces(assembly.EntrapmentSequence, digestion)
                .Where(p => p.Length >= digestion.MinLength)
                .Select(p => p.Count(c => c is 'S' or 'T')));
        }

        Assert.That(excised, Is.GreaterThan(0),
            "this fixture must actually exercise excision, or the assertion below proves nothing");
        Assert.That(entrapmentSites.OrderBy(n => n).ToList(), Is.EqualTo(expectedSites.OrderBy(n => n).ToList()),
            "apart from peptides with no possible partner, the entrapment database must present the "
            + "same distribution of candidate sites as the target");
    }

    [Test]
    public void Assemble_GivesEachFoldADifferentProtein()
    {
        const string target = "SYKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK";

        var sequences = Enumerable.Range(0, 5)
            .Select(k => EntrapmentAssembler.Assemble(target, Tryptic(), NothingForbidden, fold: k, foldCount: 5)
                .EntrapmentSequence)
            .ToList();

        Assert.That(sequences.Distinct().Count(), Is.EqualTo(5));
        Assert.That(sequences.All(s => Sorted(s) == Sorted(target)), Is.True,
            "every fold stays isomeric with the target");
    }

    // The first piece must be long enough that anchoring the protein's first two positions still
    // leaves a permutation space -- "SYK" has none once positions 0, 1 and the cleavage K are all
    // held, and the assembler then keeps it verbatim instead of permuting it.
    private const string FourPieces = "MSTQAEVDLNSGWKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK";

    /// <summary>
    /// The concatenation of the first two entrapment pieces, built with nothing forbidden. Deriving
    /// the fixture this way rather than hand-writing a sequence guarantees the collision case is
    /// actually reached: whatever the generator produces, that exact run is what gets forbidden.
    /// </summary>
    private static (EntrapmentAssembly Free, string Run) FirstRun(IDigestionParams digestion)
    {
        EntrapmentAssembly free = EntrapmentAssembler.Assemble(FourPieces, digestion, NothingForbidden);
        Assert.That(free.Pieces.Count, Is.GreaterThan(2), "fixture must have several base pieces");
        Assert.That(free.Pieces[0].Outcome, Is.EqualTo(PieceOutcome.Permuted));
        Assert.That(free.Pieces[1].Outcome, Is.EqualTo(PieceOutcome.Permuted));
        return (free, free.Pieces[0].EntrapmentPiece_ + free.Pieces[1].EntrapmentPiece_);
    }

    [Test]
    public void ARunThatWouldEqualARealPeptideIsRejected()
    {
        IDigestionParams digestion = Tryptic();
        (EntrapmentAssembly free, string run) = FirstRun(digestion);

        // Forbid exactly the missed-cleavage peptide the free build would have produced. Neither
        // piece is forbidden on its own, so only a run-aware check can avoid it.
        var forbidden = new HashSet<string> { run };
        EntrapmentAssembly guarded = EntrapmentAssembler.Assemble(FourPieces, digestion, forbidden);

        string guardedRun = guarded.Pieces[0].EntrapmentPiece_ + guarded.Pieces[1].EntrapmentPiece_;
        Assert.That(guardedRun, Is.Not.EqualTo(run), "the forbidden run must not be emitted");
        Assert.That(guarded.EntrapmentSequence, Does.Not.Contain(run));

        // The first piece has nothing before it, so no run ends at it and it is free to stay put.
        Assert.That(guarded.Pieces[0].EntrapmentPiece_, Is.EqualTo(free.Pieces[0].EntrapmentPiece_));

        // Avoiding the collision must not cost isomerism: the replacement is still a rearrangement.
        Assert.That(Sorted(guarded.Pieces[1].EntrapmentPiece_),
            Is.EqualTo(Sorted(free.Pieces[1].EntrapmentPiece_)));
        Assert.That(Sorted(guardedRun), Is.EqualTo(Sorted(run)));
    }

    [Test]
    public void RunCheckingIsOffWhenTheDigestionAllowsNoMissedCleavages()
    {
        // With no missed cleavages there are no runs to collide, so forbidding a concatenation must
        // change nothing -- otherwise the check is firing where it has no business.
        IDigestionParams noMissed = new DigestionParams("trypsin", minPeptideLength: 7, maxMissedCleavages: 0);
        EntrapmentAssembly free = EntrapmentAssembler.Assemble(FourPieces, noMissed, NothingForbidden);
        string run = free.Pieces[0].EntrapmentPiece_ + free.Pieces[1].EntrapmentPiece_;

        EntrapmentAssembly guarded = EntrapmentAssembler.Assemble(FourPieces, noMissed,
            new HashSet<string> { run });

        Assert.That(guarded.EntrapmentSequence, Is.EqualTo(free.EntrapmentSequence));
    }

    [Test]
    public void RunCheckingKeepsTheConstructionDeterministic()
    {
        IDigestionParams digestion = Tryptic();
        (_, string run) = FirstRun(digestion);
        var forbidden = new HashSet<string> { run };

        EntrapmentAssembly first = EntrapmentAssembler.Assemble(FourPieces, digestion, forbidden);
        EntrapmentAssembly second = EntrapmentAssembler.Assemble(FourPieces, digestion, forbidden);

        Assert.That(second.EntrapmentSequence, Is.EqualTo(first.EntrapmentSequence));
    }

    [Test]
    public void ACollisionOnAPieceWithNoAlternativeIsCountedNotRepaired()
    {
        // "AAR" holds its R as a cleavage site and has two identical A's left, so it has exactly one
        // arrangement. When a run ending there equals a real peptide there is no candidate to move
        // to, and the project's rule is to count such a collision rather than repair it by
        // backtracking or excising a good piece.
        IDigestionParams digestion = Tryptic();
        const string target = "MSTQAEVDLNSGWKAAR";

        EntrapmentAssembly free = EntrapmentAssembler.Assemble(target, digestion, NothingForbidden);
        Assert.That(free.Pieces[1].Outcome, Is.EqualTo(PieceOutcome.KeptVerbatimTooShort),
            "fixture must end in a piece with no alternative arrangement");
        Assert.That(free.UnrepairableRunCollisions, Is.Zero);

        string run = free.Pieces[0].EntrapmentPiece_ + free.Pieces[1].EntrapmentPiece_;
        EntrapmentAssembly guarded = EntrapmentAssembler.Assemble(target, digestion,
            new HashSet<string> { run });

        Assert.That(guarded.UnrepairableRunCollisions, Is.EqualTo(1));
        Assert.That(guarded.EntrapmentSequence, Is.EqualTo(free.EntrapmentSequence),
            "nothing is silently repaired -- the sequence is unchanged and the collision is reported");
    }

    [Test]
    public void ARunSpanningThreePiecesIsCheckedToo()
    {
        // maxMissedCleavages = 2 means a peptide can span three base pieces, so the check must look
        // back two, not one. Forbidding the three-piece run exercises the deeper lookback.
        IDigestionParams digestion = Tryptic();
        EntrapmentAssembly free = EntrapmentAssembler.Assemble(FourPieces, digestion, NothingForbidden);
        string longRun = free.Pieces[0].EntrapmentPiece_ + free.Pieces[1].EntrapmentPiece_
                         + free.Pieces[2].EntrapmentPiece_;

        EntrapmentAssembly guarded = EntrapmentAssembler.Assemble(FourPieces, digestion,
            new HashSet<string> { longRun });

        Assert.That(guarded.EntrapmentSequence, Does.Not.Contain(longRun));
        Assert.That(Sorted(guarded.Pieces[2].EntrapmentPiece_),
            Is.EqualTo(Sorted(free.Pieces[2].EntrapmentPiece_)));
    }
}
