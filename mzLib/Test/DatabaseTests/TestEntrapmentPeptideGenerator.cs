using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Numerics;
using NUnit.Framework;
using Omics.Digestion;
using UsefulProteomicsDatabases;
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

    [Test]
    public void ANullForbiddenSetMeansNothingIsForbidden()
    {
        // The parameter is optional in spirit -- a caller generating against no target database has
        // nothing to forbid -- and a null set used to reach `.Contains` as a NullReferenceException
        // out of a public API, which names no argument.
        const string target = "ALADQMNLLLSK";

        EntrapmentPeptide viaNull = EntrapmentPeptideGenerator.Create(target, Trypsin, null);
        EntrapmentPeptide viaEmpty = EntrapmentPeptideGenerator.Create(target, Trypsin,
            new HashSet<string>());

        Assert.That(viaNull.Succeeded, Is.True);
        Assert.That(viaNull.EntrapmentSequence, Is.EqualTo(viaEmpty.EntrapmentSequence),
            "null and empty must be the same request, not merely both non-throwing");
    }

    [Test]
    public void AFoldCountTheSpaceCannotCoverIsNeverReportedAsACrowdedDatabase()
    {
        // The identity is never a usable partner, so the space to share out is `size - 1`, not
        // `size`. Testing `size / foldCount == 0` missed the boundary: whichever fold's stretch
        // held the identity had its one candidate refused and reported as AllPermutationsTaken --
        // which sends a caller after a different TARGET DATABASE when the answer is a smaller fold
        // count. With nothing forbidden, AllPermutationsTaken is never the honest answer: nothing
        // took them.
        const string target = "AAGK";   // AAG, AGA, GAA with the K pinned

        Assert.That(UsefulProteomicsDatabases.DecoySequenceValidator.PermutationSpaceSize(
                target, Trypsin, null),
            Is.EqualTo(new BigInteger(3)),
            "fixture must really have a space of three, or it proves nothing");

        Assert.That(EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden).Succeeded,
            Is.True, "one fold is served by either of the two non-identity arrangements");

        for (int foldCount = 2; foldCount <= 4; foldCount++)
        {
            for (int fold = 0; fold < foldCount; fold++)
            {
                EntrapmentPeptide result = EntrapmentPeptideGenerator.Create(target, Trypsin,
                    NothingForbidden, fold: fold, foldCount: foldCount);

                Assert.That(result.Failure, Is.Not.EqualTo(EntrapmentFailure.AllPermutationsTaken),
                    "fold " + fold + " of " + foldCount + ": nothing was forbidden, so no "
                    + "arrangement can have been taken -- the fold count is what does not fit");

                if (!result.Succeeded)
                {
                    Assert.That(result.Failure,
                        Is.EqualTo(EntrapmentFailure.SpaceTooSmallForFoldCount));
                }
            }
        }

        // Three folds over three arrangements cannot work at all, one of the three being the target.
        for (int fold = 0; fold < 3; fold++)
        {
            Assert.That(EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden,
                    fold: fold, foldCount: 3).Failure,
                Is.EqualTo(EntrapmentFailure.SpaceTooSmallForFoldCount));
        }
    }

    [Test]
    public void AFailedSearchReportsHowMuchOfTheSpaceItWalked()
    {
        // ProbesUsed was hard-coded to 0 on every failure path, so the one number saying the walk
        // was EXHAUSTIVE -- which is what makes "no partner exists" a proof rather than an
        // abandoned search -- was absent exactly where it mattered.
        const string target = "AAGK";
        var everythingTaken = new HashSet<string> { "AGAK", "GAAK" };

        EntrapmentPeptide result = EntrapmentPeptideGenerator.Create(target, Trypsin, everythingTaken);

        Assert.That(result.Succeeded, Is.False);
        Assert.That(result.Failure, Is.EqualTo(EntrapmentFailure.AllPermutationsTaken),
            "here the database really did take them, which is the other side of the test above");
        Assert.That(result.ProbesUsed, Is.EqualTo(2),
            "both non-identity arrangements were examined, and the count is the evidence of it; "
            + "the identity is never in a fold's share, so it is not probed");
    }

    [Test]
    public void EveryFoldGetsACandidateWheneverTheSpaceCoversTheFoldCount()
    {
        // Alexander-Sol's repro on #1271. With the termini anchored, MAAABAK frees AAABA: five
        // arrangements, the identity at rank 1. Contiguous stretches of size / foldCount = 1 gave
        // fold 1 only the identity and left ranks 3 and 4 unused, so fold 1 was excised as
        // SpaceTooSmallForFoldCount although four non-identity arrangements existed for three folds.
        const string target = "MAAABAK";
        int[] anchors = { 0, target.Length - 1 };

        Assert.That(DecoySequenceValidator.PermutationSpaceSize(target, Trypsin, anchors),
            Is.EqualTo(new BigInteger(5)), "fixture must really have a space of five");
        Assert.That(DecoySequenceValidator.RankPermutation(target, Trypsin, anchors),
            Is.EqualTo(BigInteger.One), "fixture must put the identity at rank 1");

        var partners = new HashSet<string>();
        for (int seed = 1; seed <= 20; seed++)
        {
            partners.Clear();
            for (int fold = 0; fold < 3; fold++)
            {
                EntrapmentPeptide result = EntrapmentPeptideGenerator.Create(target, Trypsin,
                    NothingForbidden, fold: fold, foldCount: 3, seed: seed, alsoHeldInPlace: anchors);

                Assert.That(result.Succeeded, Is.True, $"seed {seed}, fold {fold}: {result.Failure}");
                Assert.That(result.EntrapmentSequence, Is.Not.EqualTo(target));
                partners.Add(result.EntrapmentSequence!);
            }

            Assert.That(partners, Has.Count.EqualTo(3), "folds draw on disjoint shares, so never repeat");
        }

        // And the boundary from the test above: two folds over AAGK's two non-identity arrangements.
        for (int fold = 0; fold < 2; fold++)
        {
            Assert.That(EntrapmentPeptideGenerator.Create("AAGK", Trypsin, NothingForbidden,
                fold: fold, foldCount: 2).Succeeded, Is.True);
        }
    }

    [Test]
    public void ASingleFoldIsNotConfinedToOneBlockOfTheLexicographicOrder()
    {
        // Alexander-Sol, #1271: with the termini anchored, AEGLSVTK frees EGLSVT -- 720
        // arrangements. Contiguous stretches gave fold 0 of r = 9 only the first 80, so it always
        // put E second, and fold 8 always put V there, whatever the seed. A per-fold estimate then
        // sampled a skewed population. Interleaved shares span the whole order, so the seed decides.
        const string target = "AEGLSVTK";
        int[] anchors = { 0, target.Length - 1 };
        for (int fold = 0; fold < 9; fold++)
        {
            var secondResidues = new HashSet<char>();
            for (int seed = 1; seed <= 40; seed++)
            {
                EntrapmentPeptide result = EntrapmentPeptideGenerator.Create(target, Trypsin,
                    NothingForbidden, fold: fold, foldCount: 9, seed: seed, alsoHeldInPlace: anchors);
                Assert.That(result.Succeeded, Is.True);
                secondResidues.Add(result.EntrapmentSequence![1]);
            }

            Assert.That(secondResidues.Count, Is.GreaterThan(1),
                $"fold {fold} put the same residue second under every seed");
        }
    }

    [Test]
    public void ASingleFoldIsNotConfinedToOneOrderOfItsLastFreeResidues()
    {
        // Alexander-Sol, #1271 third pass: interleaving by residue class modulo r moved the per-fold
        // skew to the other end. The low digits of a lexicographic rank encode the relative order of
        // the LAST free positions, so with nonIdentityRank = fold + r * j, rank mod 2 fixed whether the
        // last two were ascending and rank mod 6 the order of the last three: at r = 6 every fold of
        // AEGLSVTK saw one order of positions 4-6 under all 200 seeds, next to the pinned K. A keyed
        // bijection over the non-identity ranks, applied before the residue-class split, breaks it.
        const string target = "AEGLSVTK";
        int[] anchors = { 0, target.Length - 1 };
        foreach (int foldCount in new[] { 2, 6, 9 })
        {
            for (int fold = 0; fold < foldCount; fold++)
            {
                var tailOrders = new HashSet<string>();
                var lastTwoAscending = new HashSet<bool>();
                for (int seed = 1; seed <= 200; seed++)
                {
                    EntrapmentPeptide result = EntrapmentPeptideGenerator.Create(target, Trypsin,
                        NothingForbidden, fold: fold, foldCount: foldCount, seed: seed, alsoHeldInPlace: anchors);
                    Assert.That(result.Succeeded, Is.True);
                    string tail = result.EntrapmentSequence!.Substring(4, 3);
                    tailOrders.Add(string.Concat(tail.Select(c => tail.Count(o => o < c))));
                    lastTwoAscending.Add(tail[1] < tail[2]);
                }

                Assert.That(tailOrders.Count, Is.GreaterThan(1),
                    $"r = {foldCount}, fold {fold} put its last three free residues in one order under every seed");
                Assert.That(lastTwoAscending, Has.Count.EqualTo(2),
                    $"r = {foldCount}, fold {fold} put its last two free residues in one order under every seed");
            }
        }
    }

    [Test]
    public void TheShuffledSharesStillPartitionEveryNonIdentityArrangement()
    {
        // The keyed shuffle has to be a bijection on the non-identity ranks, or folds would overlap
        // (no longer distinct partners) or leave arrangements no fold can reach (a false proof of
        // exhaustion). AEGLSK anchored frees EGLS: 24 arrangements,
        // 23 usable, so the shuffle works over 5 bits and must cycle-walk 32 back into 23 with an
        // uneven split. Forbid every arrangement but one: exactly one fold must find it.
        const string target = "AEGLSK";
        int[] anchors = { 0, target.Length - 1 };
        BigInteger size = DecoySequenceValidator.PermutationSpaceSize(target, Trypsin, anchors);
        Assert.That(size, Is.EqualTo(new BigInteger(24)), "fixture must really have a space of 24");

        var arrangements = new List<string>();
        for (BigInteger i = BigInteger.Zero; i < size; i++)
        {
            string arrangement = DecoySequenceValidator.UnrankPermutation(target, Trypsin, i, out _, anchors);
            if (arrangement != target)
            {
                arrangements.Add(arrangement);
            }
        }

        foreach (int foldCount in new[] { 1, 3 })
        {
            for (int seed = 1; seed <= 3; seed++)
            {
                foreach (string onlyFree in arrangements)
                {
                    var forbidden = new HashSet<string>(arrangements.Where(a => a != onlyFree));
                    int foundBy = Enumerable.Range(0, foldCount).Count(fold => EntrapmentPeptideGenerator.Create(
                        target, Trypsin, forbidden, fold: fold, foldCount: foldCount, seed: seed,
                        alsoHeldInPlace: anchors).EntrapmentSequence == onlyFree);

                    Assert.That(foundBy, Is.EqualTo(1), $"r = {foldCount}, seed {seed}: {onlyFree}");
                }
            }
        }
    }

    [Test]
    public void RankPermutationIsTheInverseOfUnrankAtTheIdentity()
    {
        string[] sequences = { "AEGLSVTK", "MAAABAK", "AAGK", "PEPTIDEKAAR", "SSSSSSR", "K", "LIHTGVKLIHTVGK" };
        foreach (string sequence in sequences)
        {
            foreach (int[]? anchors in new[] { null, new[] { 0, sequence.Length - 1 } })
            {
                BigInteger rank = DecoySequenceValidator.RankPermutation(sequence, Trypsin, anchors);
                string back = DecoySequenceValidator.UnrankPermutation(sequence, Trypsin, rank, out _, anchors);
                Assert.That(back, Is.EqualTo(sequence), sequence);
            }
        }
    }

    // ---- production review of #1271, 2026-10-10 --------------------------------

    [TestCase("K[P]|,R[P]|", "preventing-cleavage motif")]
    [TestCase("TX|T", "multi-residue motif")]
    [TestCase("X|", "wildcard motif")]
    public void Create_RefusesMotifsWhoseSitesItCannotHold(string motifs, string named)
    {
        // The refusal lived only on the protein path, so this public method handed back partners
        // that digest differently from their targets.
        var thrown = Assert.Throws<MzLibUtil.MzLibException>(() => EntrapmentPeptideGenerator.Create(
            "AKLPPR", DigestionMotif.ParseDigestionMotifsFromString(motifs), NothingForbidden));
        Assert.That(thrown!.Message, Does.Contain(named));
    }

    [Test]
    public void Create_PastTheRefusalAPreventingMotifReallyMovesACleavageSite()
    {
        // Why the refusal is not paranoia. Under trypsin|P neither the K nor the P after it is held,
        // so the unvetted path is free to move the P next to the K and erase the site, or away from
        // a K it was protecting and invent one.
        List<DigestionMotif> trypsinP = DigestionMotif.ParseDigestionMotifsFromString("K[P]|,R[P]|");
        HashSet<int> sites = DecoySequenceValidator.CleavageSitePositions("AKLPPR", trypsinP);

        bool aSiteMoved = Enumerable.Range(0, 2).Any(fold =>
        {
            EntrapmentPeptide partner = EntrapmentPeptideGenerator.CreateFromVettedMotifs("AKLPPR", trypsinP,
                NothingForbidden, fold, foldCount: 2, seed: 1, alsoHeldInPlace: new[] { 0, 5 },
                rejectInContext: null);
            return !DecoySequenceValidator.CleavageSitePositions(partner.EntrapmentSequence!, trypsinP)
                .SetEquals(sites);
        });

        Assert.That(aSiteMoved, Is.True);
    }

    [Test]
    public void Create_RefusesMissingOrEmptyMotifs()
    {
        var empty = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentPeptideGenerator.Create("AKLPPR", new List<DigestionMotif>(), NothingForbidden));
        Assert.That(empty!.Message, Does.Contain("no cleavage motifs"));

        var missing = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentPeptideGenerator.Create("AKLPPR", null!, NothingForbidden));
        Assert.That(missing!.Message, Does.Contain("cleavage motifs"));
    }

    [Test]
    public void Create_ReportsARefusalByContextAsARunCollisionNotAsATakenSpace()
    {
        // The two failures send a caller after different remedies -- a different seed or different
        // neighbours, against a different target database -- and no test checked that a refusal by
        // context was reported as the first rather than the second.
        EntrapmentPeptide result = EntrapmentPeptideGenerator.Create("ACDEFGHIK", Trypsin, NothingForbidden,
            rejectInContext: _ => true);

        Assert.That(result.Succeeded, Is.False);
        Assert.That(result.Failure, Is.EqualTo(EntrapmentFailure.RunCollisionsExhaustedTheSpace));
        Assert.That(result.ProbesUsed, Is.EqualTo((int)(result.PermutationSpaceSize - 1)),
            "the whole share was walked, which is what makes the failure a proof");
    }

    [Test]
    public void Create_EachFoldReachesItsWholeShareEvenlyAcrossSeeds()
    {
        // ASingleFoldIsNotConfinedToOneOrderOfItsLastFreeResidues asks only for "more than one order"
        // per fold, which a 199:1 skew passes. This asks for the distribution. AEGLSK with its ends
        // held has 23 non-identity arrangements, so over 2,300 seeds every fold should reach every
        // one about 100 times; residue classes without the keyed shuffle gave each fold a third.
        const string target = "AEGLSK";
        int[] anchors = { 0, 5 };
        for (int fold = 0; fold < 3; fold++)
        {
            var counts = new Dictionary<string, int>();
            for (int seed = 1; seed <= 2300; seed++)
            {
                string partner = EntrapmentPeptideGenerator.Create(target, Trypsin, NothingForbidden,
                    fold: fold, foldCount: 3, seed: seed, alsoHeldInPlace: anchors).EntrapmentSequence!;
                counts[partner] = counts.GetValueOrDefault(partner) + 1;
            }

            Assert.That(counts.Count, Is.EqualTo(23), $"fold {fold} reaches every non-identity arrangement");
            Assert.That(counts.Values.Min(), Is.GreaterThan(50), $"fold {fold}: no arrangement starved");
            Assert.That(counts.Values.Max(), Is.LessThan(150), $"fold {fold}: no arrangement favoured");
        }
    }

    [Test]
    public void Create_PartnersArePinnedToLiteralSequences()
    {
        // Every other determinism test compares the generator with itself in one process, so a change
        // to which partner is chosen -- the key, the round function, the offset, the fold allocation
        // -- passed them all while every database anyone had built stopped regenerating. These
        // literals are construction 1's. If this fails, partners have changed: bump
        // EntrapmentProteinGenerator.ConstructionVersion and say so in the release notes, then update
        // the literals. The negative seed covers the invariant formatting of the key material.
        Assert.That(EntrapmentProteinGenerator.ConstructionVersion, Is.EqualTo(1),
            "the literals below belong to construction 1");

        Assert.That(EntrapmentPeptideGenerator.Create("SYKALADQMNLLLSK", Trypsin, NothingForbidden, seed: 1)
            .EntrapmentSequence, Is.EqualTo("SAKLMQLLSLYNDAK"));

        string[] folds = Enumerable.Range(0, 3).Select(fold => EntrapmentPeptideGenerator.Create(
                "SYKALADQMNLLLSK", Trypsin, NothingForbidden, fold: fold, foldCount: 3, seed: -7,
                alsoHeldInPlace: new[] { 0, 14 }).EntrapmentSequence!)
            .ToArray();
        Assert.That(folds, Is.EqualTo(new[] { "SAKNQLLDMLAYLSK", "SAKMYLDLAQLSNLK", "SNKAASLLLMLYQDK" }));
    }
}
