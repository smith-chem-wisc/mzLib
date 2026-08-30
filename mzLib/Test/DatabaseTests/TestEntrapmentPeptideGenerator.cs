using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Numerics;
using NUnit.Framework;
using Omics.Digestion;
using UsefulProteomicsDatabases.EntrapmentGeneration;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class EntrapmentPeptideGeneratorTests
{
    private static List<DigestionMotif> Trypsin => DigestionMotif.ParseDigestionMotifsFromString("K|,R|");
    private static readonly HashSet<string> NothingForbidden = new();

    [Test]
    public void Create_ReturnsAnIsomericPeptideThatIsNotTheTarget()
    {
        EntrapmentPeptide result = EntrapmentPeptideGenerator.Create("SYKALADQMNLLLSK", Trypsin, NothingForbidden);

        Assert.That(result.Succeeded, Is.True);
        Assert.That(result.Failure, Is.EqualTo(EntrapmentFailure.None));
        Assert.That(result.EntrapmentSequence, Is.Not.EqualTo("SYKALADQMNLLLSK"));
        Assert.That(string.Concat(result.EntrapmentSequence!.OrderBy(c => c)),
            Is.EqualTo(string.Concat("SYKALADQMNLLLSK".OrderBy(c => c))),
            "same residues in a different order -- same mass, same composition");
    }

    [Test]
    public void Create_PreservesCandidateSiteCountAndCleavageSites()
    {
        // The property the glyco localization work depends on: the number of S/T in a peptide is
        // the number of candidate sites, and comparing site calls at different n compares
        // different difficulty. Composition preservation gives this for free -- assert it anyway.
        const string target = "SYKALADQMNLLLSK";
        EntrapmentPeptide result = EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden);

        Assert.That(result.EntrapmentSequence!.Count(c => c is 'S' or 'T'),
            Is.EqualTo(target.Count(c => c is 'S' or 'T')));
        Assert.That(result.EntrapmentSequence[2], Is.EqualTo('K'));
        Assert.That(result.EntrapmentSequence[14], Is.EqualTo('K'));
    }

    [Test]
    public void Create_IsAPureFunctionOfSequenceAndSeed()
    {
        // No random state, so the answer cannot depend on how many peptides were processed first,
        // on which protein this one came from, or on which thread ran it.
        const string target = "SYKALADQMNLLLSK";
        string first = EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden).EntrapmentSequence!;

        for (int i = 0; i < 5; i++)
        {
            EntrapmentPeptideGenerator.Create("ACDEFGHIK", Trypsin, NothingForbidden);   // unrelated work
        }

        Assert.That(EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden).EntrapmentSequence,
            Is.EqualTo(first));
        Assert.That(EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden, seed: 2).EntrapmentSequence,
            Is.Not.EqualTo(first), "a different seed must give a different answer");
    }

    [Test]
    public void Create_ReportsNoPermutationExistsSeparatelyFromExhaustion()
    {
        // Two different failures that must never be merged into one "failed":
        //   arithmetic -- a homopolymeric tract has exactly one arrangement, so no r helps;
        //   exhaustion -- the space is real but every member of it is already spoken for.
        EntrapmentPeptide homopolymer = EntrapmentPeptideGenerator.Create("SSSSSSR", Trypsin, NothingForbidden);
        Assert.That(homopolymer.Succeeded, Is.False);
        Assert.That(homopolymer.Failure, Is.EqualTo(EntrapmentFailure.NoPermutationExists));
        Assert.That(homopolymer.PermutationSpaceSize, Is.EqualTo(BigInteger.One));

        // Forbid the entire space of a small peptide, leaving exhaustion as the only outcome.
        var everyPermutation = new HashSet<string>();
        BigInteger size = UsefulProteomicsDatabases.DecoySequenceValidator.PermutationSpaceSize("QEEEKKK", Trypsin);
        for (BigInteger i = BigInteger.Zero; i < size; i++)
        {
            everyPermutation.Add(UsefulProteomicsDatabases.DecoySequenceValidator
                .UnrankPermutation("QEEEKKK", Trypsin, i, out _));
        }

        EntrapmentPeptide exhausted = EntrapmentPeptideGenerator.Create("QEEEKKK", Trypsin, everyPermutation);
        Assert.That(exhausted.Succeeded, Is.False);
        Assert.That(exhausted.Failure, Is.EqualTo(EntrapmentFailure.AllPermutationsTaken));
        Assert.That(exhausted.PermutationSpaceSize, Is.GreaterThan(BigInteger.One),
            "the space exists -- it is simply fully occupied");
    }

    [Test]
    public void Create_ReportsWhenTheSpaceCannotSupplyTheRequestedFolds()
    {
        // "QEEEKKK" has four arrangements. Asking for nine folds is not exhaustion and not
        // arithmetic impossibility -- it is a request the database cannot honour, and the caller's
        // remedy (ask for fewer folds) differs from both.
        EntrapmentPeptide result = EntrapmentPeptideGenerator.Create("QEEEKKK", Trypsin, NothingForbidden, fold: 8, foldCount: 9);

        Assert.That(result.Succeeded, Is.False);
        Assert.That(result.Failure, Is.EqualTo(EntrapmentFailure.SpaceTooSmallForFoldCount));
    }

    [Test]
    public void Create_GivesMutuallyDistinctPeptidesForEachFold()
    {
        const string target = "SYKALADQMNLLLSK";
        const int folds = 9;

        var sequences = Enumerable.Range(0, folds)
            .Select(k => EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden, fold: k, foldCount: folds))
            .ToList();

        Assert.That(sequences.All(s => s.Succeeded), Is.True);
        Assert.That(sequences.Select(s => s.EntrapmentSequence).Distinct().Count(), Is.EqualTo(folds),
            "r-fold entrapment needs r distinct partners, not r copies");
        Assert.That(sequences.All(s => s.EntrapmentSequence != target), Is.True);
    }

    [Test]
    public void Create_FoldsAreIndependentOfOneAnother()
    {
        // Fold 5 must be the same whether or not folds 0-4 were ever asked for, so a database can
        // be extended or regenerated in pieces without invalidating what already exists.
        const string target = "SYKALADQMNLLLSK";

        string alone = EntrapmentPeptideGenerator
            .Create(target, Trypsin, NothingForbidden, fold: 5, foldCount: 9).EntrapmentSequence!;

        for (int k = 0; k < 5; k++)
        {
            EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden, fold: k, foldCount: 9);
        }

        Assert.That(EntrapmentPeptideGenerator
            .Create(target, Trypsin, NothingForbidden, fold: 5, foldCount: 9).EntrapmentSequence,
            Is.EqualTo(alone));
    }

    [Test]
    public void Create_NeverReturnsAForbiddenSequence()
    {
        const string target = "ACDEFGHIK";
        var forbidden = new HashSet<string>();

        // Forbid the first few answers and confirm the generator walks past them.
        for (int i = 0; i < 3; i++)
        {
            EntrapmentPeptide step = EntrapmentPeptideGenerator.Create(target, Trypsin, forbidden);
            Assert.That(step.Succeeded, Is.True);
            Assert.That(forbidden, Does.Not.Contain(step.EntrapmentSequence));
            forbidden.Add(step.EntrapmentSequence!);
        }

        Assert.That(forbidden.Count, Is.EqualTo(3));
    }

    [Test]
    public void Create_ReportsTheSpaceSizeAndHowHardItHadToLook()
    {
        EntrapmentPeptide easy = EntrapmentPeptideGenerator.Create("TTTPAPTTT", Trypsin, NothingForbidden);

        Assert.That(easy.PermutationSpaceSize, Is.EqualTo(new BigInteger(252)));
        Assert.That(easy.ProbesUsed, Is.EqualTo(1), "an unforbidden space should answer on the first probe");
    }

    // Each guard is asserted through its MESSAGE, not merely through the exception type. Several of
    // these arguments are rejected by more than one guard, so a type-only assertion still passes
    // when the guard under test is removed -- the message is what pins which one actually fired.
    [Test]
    public void Create_RejectsAnEmptySequence()
    {
        var thrown = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentPeptideGenerator.Create("", Trypsin, NothingForbidden));
        Assert.That(thrown!.Message, Does.Contain("empty"));
    }

    [Test]
    public void Create_RejectsAFoldCountBelowOne()
    {
        var thrown = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentPeptideGenerator.Create("ACDEFGHIK", Trypsin, NothingForbidden, fold: 0, foldCount: 0));
        Assert.That(thrown!.Message, Does.Contain("Fold count").And.Contain("0"));
    }

    [Test]
    public void Create_RejectsAFoldOutsideTheRequestedCount()
    {
        var tooHigh = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentPeptideGenerator.Create("ACDEFGHIK", Trypsin, NothingForbidden, fold: 3, foldCount: 3));
        Assert.That(tooHigh!.Message, Does.Contain("Fold 3").And.Contain("3 requested folds"));

        var negative = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentPeptideGenerator.Create("ACDEFGHIK", Trypsin, NothingForbidden, fold: -1));
        Assert.That(negative!.Message, Does.Contain("Fold -1"));
    }
}
