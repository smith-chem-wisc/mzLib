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
        Assert.That(assembly.KeptVerbatimCount, Is.GreaterThan(0));
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

        Assert.That(assembly.ExcisedCount, Is.EqualTo(1));
        Assert.That(assembly.MissedCleavagePeptidesSpanningAnExcision, Is.GreaterThan(0),
            "an excision in the middle of a protein necessarily breaks the peptides that spanned it");
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
}
