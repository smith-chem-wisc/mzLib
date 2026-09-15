using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// What <c>Protein.Digest</c> returns for every combination of <see cref="DigestionParams.SearchModeType"/> and
    /// <see cref="DigestionParams.FragmentationTerminus"/>. Read this first if those two settings are confusing: they are.
    /// </summary>
    /// <remarks>
    /// <para><b>The table.</b> For a fully specific protease such as trypsin:</para>
    /// <code>
    ///   SearchModeType  FragmentationTerminus  Protein.Digest returns                   answered by
    ///   Full            Both, N or C           fully specific peptides                  ProteinDigestion.Digestion
    ///   Semi            Both (the default)     semi-specific peptides                   ProteinDigestion.SemiSpecificDigestion
    ///   Semi            N                      N-terminal SEEDS, not peptides           ProteinDigestion.SpeedySemiSpecificDigestion
    ///   Semi            C                      C-terminal SEEDS, not peptides           ProteinDigestion.SpeedySemiSpecificDigestion
    ///   None            N                      singleN SEEDS (non-specific), not peptides   the singleN protease
    ///   None            C                      singleC SEEDS (non-specific), not peptides   the singleC protease
    ///   None            Both                   singleC seeds, the same as None + C      the singleC protease
    /// </code>
    ///
    /// <para><b>Peptides versus seeds.</b> A peptide is a candidate a search engine can score as it is: its mass is its
    /// mass. A seed is a long stretch fixed at ONE terminus whose other end has not been decided. Only an engine that
    /// decides that end afterwards may ask for seeds. In MetaMorpheus that is the non-specific search engine
    /// (<c>NonSpecificEnzymeSearchEngine</c>): it scores a seed with fragment ions from the fixed terminus only, then
    /// walks along the seed adding residue masses until the total matches the precursor mass, and cuts there. That works
    /// because the end point is the only unknown. It runs one pass with N and one with C, which is why N and C mean
    /// seeds.</para>
    ///
    /// <para><b>Every other engine needs peptides.</b> Classic, Modern, Glyco and crosslink searches use
    /// FragmentationTerminus Both, the default. Glyco in particular cannot use seeds: it computes precursor mass minus
    /// peptide mass to find the glycan, which is undefined when the peptide's end is not known. Before #1303, Semi + Both
    /// returned C-terminal seeds, and O-glyco searches configured that way lost most of their identifications without an
    /// error.</para>
    ///
    /// <para><b>The one sharp edge left.</b> None + Both returns singleC seeds (the constructor picks singleC for
    /// anything that is not N). No engine that needs peptides should ask for None; there is no "non-specific peptides"
    /// request. <see cref="Digest_NoneWithBothTermini_ReturnsTheSingleCSeedsNotNonSpecificPeptides"/> pins this so a
    /// change to it is deliberate.</para>
    ///
    /// <para><b>A protease whose own specificity is Semi</b> (a user-defined semi-trypsin, say) returns semi-specific
    /// peptides with SearchModeType Full; see <c>SemiSpecificDigestionTests</c>, which also checks every peptide of the
    /// Semi rows against a brute-force answer key.</para>
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class SearchModeTypeDigestionTests
    {
        /// <summary>
        /// Initiator Met, several K/R, and no K or R followed by P, so the result does not depend on whether the embedded
        /// "trypsin" applies the proline rule (its naming is being changed by #1186).
        /// </summary>
        private const string Sequence = "MAAKCCKDDKEEKFFKGGPEPTIDERSSTTAAK";

        /// <summary>What a combination of settings asks <c>Protein.Digest</c> for.</summary>
        public enum WhatDigestReturns
        {
            FullySpecificPeptides,
            SemiSpecificPeptides,
            NTerminalSeeds,
            CTerminalSeeds,
            SingleNSeeds,
            SingleCSeeds,
        }

        private static DigestionParams Params(CleavageSpecificity searchModeType, FragmentationTerminus terminus) =>
            new("trypsin", maxMissedCleavages: 2, minPeptideLength: 5, maxPeptideLength: 30,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable,
                searchModeType: searchModeType, fragmentationTerminus: terminus);

        private static List<(int Start, int End, CleavageSpecificity Label)> Digest(Protein protein, DigestionParams digestionParams) =>
            protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Select(p => (p.OneBasedStartResidue, p.OneBasedEndResidue, p.CleavageSpecificityForFdrCategory))
                .OrderBy(p => p).ToList();

        private static List<(int Start, int End, CleavageSpecificity Label)> Peptides(IEnumerable<ProteolyticPeptide> peptides) =>
            peptides.Select(p => (p.OneBasedStartResidue, p.OneBasedEndResidue, p.CleavageSpecificityForFdrCategory))
                .OrderBy(p => p).ToList();

        private static bool IsTrypticSite(int oneBasedResidue) => Sequence[oneBasedResidue - 1] is 'K' or 'R';

        /// <summary>
        /// Each row of the table in the class remarks: <c>Protein.Digest</c> gives exactly what the named
        /// <see cref="ProteinDigestion"/> method gives, and the result has the defining property of its kind.
        /// </summary>
        /// <remarks>
        /// The first check pins WHICH code answers each request, so a change in routing (the #1303 bug was one) fails
        /// here by name. The second pins what that answer MEANS, independently of the code that produced it.
        /// </remarks>
        [Test]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.Both, WhatDigestReturns.FullySpecificPeptides)]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.N, WhatDigestReturns.FullySpecificPeptides)]
        [TestCase(CleavageSpecificity.Full, FragmentationTerminus.C, WhatDigestReturns.FullySpecificPeptides)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.Both, WhatDigestReturns.SemiSpecificPeptides)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.N, WhatDigestReturns.NTerminalSeeds)]
        [TestCase(CleavageSpecificity.Semi, FragmentationTerminus.C, WhatDigestReturns.CTerminalSeeds)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.N, WhatDigestReturns.SingleNSeeds)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.C, WhatDigestReturns.SingleCSeeds)]
        [TestCase(CleavageSpecificity.None, FragmentationTerminus.Both, WhatDigestReturns.SingleCSeeds)]
        public static void Digest_EachSearchModeTypeAndTerminus_ReturnsWhatTheTableSays(CleavageSpecificity searchModeType, FragmentationTerminus terminus, WhatDigestReturns expected)
        {
            Assert.That(Sequence.Contains("KP") || Sequence.Contains("RP"), Is.False, "the sequence must not depend on trypsin's proline rule");
            var protein = new Protein(Sequence, "SearchModeTypeDigestionTests");
            DigestionParams digestionParams = Params(searchModeType, terminus);
            var digestion = new ProteinDigestion(digestionParams, new List<Modification>(), new List<Modification>());
            var actual = Digest(protein, digestionParams);

            Assert.That(actual, Is.Not.Empty, "the test protein must produce something for every row");
            Assert.That(digestionParams.SpecificProtease.Name, Is.EqualTo("trypsin"), "SpecificProtease is always the protease the caller named");

            var fullySpecific = Digest(protein, Params(CleavageSpecificity.Full, FragmentationTerminus.Both));
            var semiSpecific = Digest(protein, Params(CleavageSpecificity.Semi, FragmentationTerminus.Both));

            switch (expected)
            {
                case WhatDigestReturns.FullySpecificPeptides:
                    Assert.That(actual, Is.EqualTo(Peptides(digestion.Digestion(protein))), "answered by ProteinDigestion.Digestion");
                    Assert.That(actual, Is.EqualTo(fullySpecific), "for Full, FragmentationTerminus does not change the peptides");
                    Assert.That(actual.Select(p => p.Label), Is.All.EqualTo(CleavageSpecificity.Full));
                    break;

                case WhatDigestReturns.SemiSpecificPeptides:
                    Assert.That(actual, Is.EqualTo(Peptides(digestion.SemiSpecificDigestion(protein))), "answered by ProteinDigestion.SemiSpecificDigestion");
                    Assert.That(actual.Select(p => (p.Start, p.End)), Is.SupersetOf(fullySpecific.Select(p => (p.Start, p.End))), "semi-specific peptides include the fully specific ones");
                    Assert.That(actual.Any(p => p.Label == CleavageSpecificity.Semi), Is.True, "and add peptides with one ragged end");
                    break;

                case WhatDigestReturns.NTerminalSeeds:
                case WhatDigestReturns.CTerminalSeeds:
                    Assert.That(actual, Is.EqualTo(Peptides(digestion.SpeedySemiSpecificDigestion(protein))), "answered by ProteinDigestion.SpeedySemiSpecificDigestion");
                    Assert.That(ProteinDigestion.WantsSemiSpecificSeeds(digestionParams), Is.True);
                    Assert.That(actual.Count, Is.LessThan(semiSpecific.Count), "seeds are far fewer than the semi-specific peptides they stand for; they are not those peptides");
                    if (expected == WhatDigestReturns.NTerminalSeeds)
                        Assert.That(actual.Where(p => !(p.Start == 1 || IsTrypticSite(p.Start - 1) || (p.Start == 2 && Sequence[0] == 'M'))), Is.Empty,
                            "every N seed starts on a specific N-terminus; its C-terminal end is left for the search engine to decide");
                    else
                        Assert.That(actual.Where(p => !(p.End == Sequence.Length || IsTrypticSite(p.End))), Is.Empty,
                            "every C seed ends on a specific C-terminus; its N-terminal end is left for the search engine to decide");
                    break;

                case WhatDigestReturns.SingleNSeeds:
                case WhatDigestReturns.SingleCSeeds:
                    string swappedProtease = expected == WhatDigestReturns.SingleNSeeds ? "singleN" : "singleC";
                    CleavageSpecificity label = expected == WhatDigestReturns.SingleNSeeds ? CleavageSpecificity.SingleN : CleavageSpecificity.SingleC;
                    Assert.That(digestionParams.Protease.Name, Is.EqualTo(swappedProtease), "SearchModeType None swaps Protease for singleN or singleC and keeps the named one as SpecificProtease");
                    Assert.That(actual, Is.EqualTo(Peptides(digestion.Digestion(protein))), $"answered by the {swappedProtease} protease's digestion");
                    Assert.That(actual.Select(p => p.Label), Is.All.EqualTo(label));
                    break;
            }
        }

        /// <summary>
        /// None + Both is not a request for non-specific peptides. It returns exactly the singleC seeds of None + C, which
        /// only an engine that trims the N-terminal end afterwards can use.
        /// </summary>
        /// <remarks>
        /// This pins current behaviour so that changing it is a deliberate decision, not an accident. A caller that needs
        /// peptides (Classic, Modern, Glyco, crosslink) must not use SearchModeType None; MetaMorpheus uses None only in
        /// its non-specific search, which always runs separate N and C passes.
        /// </remarks>
        [Test]
        public static void Digest_NoneWithBothTermini_ReturnsTheSingleCSeedsNotNonSpecificPeptides()
        {
            var protein = new Protein(Sequence, "SearchModeTypeDigestionTests");
            var noneBoth = Digest(protein, Params(CleavageSpecificity.None, FragmentationTerminus.Both));
            var noneC = Digest(protein, Params(CleavageSpecificity.None, FragmentationTerminus.C));

            Assert.That(noneBoth, Is.EqualTo(noneC));
            Assert.That(noneBoth.Select(p => p.End).Distinct().Count(), Is.EqualTo(noneBoth.Count),
                "singleC seeds: one per C-terminal end, so these cannot be the non-specific peptides (many per end)");
        }

        /// <summary>
        /// Cloning to a new terminus, which is how MetaMorpheus's non-specific search makes its N and C passes, changes
        /// the answer from peptides to seeds. Pins that a clone is a different request, not the same one.
        /// </summary>
        [Test]
        [TestCase(CleavageSpecificity.Semi)]
        [TestCase(CleavageSpecificity.None)]
        public static void Clone_ToNOrC_TurnsAPeptideRequestIntoASeedRequest(CleavageSpecificity searchModeType)
        {
            var protein = new Protein(Sequence, "SearchModeTypeDigestionTests");
            DigestionParams both = Params(searchModeType, FragmentationTerminus.Both);
            var cloneN = (DigestionParams)both.Clone(FragmentationTerminus.N);
            var cloneC = (DigestionParams)both.Clone(FragmentationTerminus.C);

            Assert.That(cloneN.SpecificProtease.Name, Is.EqualTo("trypsin"));
            Assert.That(cloneC.SpecificProtease.Name, Is.EqualTo("trypsin"));
            Assert.That(Digest(protein, cloneN), Is.EqualTo(Digest(protein, Params(searchModeType, FragmentationTerminus.N))));
            Assert.That(Digest(protein, cloneC), Is.EqualTo(Digest(protein, Params(searchModeType, FragmentationTerminus.C))));
            if (searchModeType == CleavageSpecificity.Semi)
                Assert.That(Digest(protein, cloneN), Is.Not.EqualTo(Digest(protein, both)), "seeds are not the semi-specific peptides");
        }
    }
}
