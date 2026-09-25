using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// End-to-end digestion tests for a protease whose motif REQUIRES a modification at one of its
    /// subsites. Slice 1 classified the chemistry and slice 2 acts on it: these are the first tests in
    /// which a glycan actually decides whether a peptide bond is severed.
    /// </summary>
    /// <remarks>
    /// <para>The substrates are the published ones from the truth set
    /// (<c>TestData/digestion-truth-set.tsv</c>), so the expectations come from the enzymology rather
    /// than from the implementation. StcE on <c>RPPIT*QSSL</c> is the canonical case: Malaker et al.
    /// report it "converting RPPIT*QSSL to RPPIT*Q", and report that the same backbone carrying
    /// beta-O-GlcNAc instead of alpha-O-GalNAc is not cleaved at all.</para>
    ///
    /// <para>No protease in <c>proteases.tsv</c> declares a requirement yet -- that is slice 3 -- so
    /// these fixtures build one in code. That is deliberate: it keeps the mechanism under test separate
    /// from the data file that will eventually configure it, and it lets the requirement be varied
    /// (P2 vs P1', class, subsite number) in ways a shipped entry never would be.</para>
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CleavagePromotingDigestionTests
    {
        /// <summary>An O-linked glycan on Ser/Thr, classified through the ModificationType MetaMorpheus sets.</summary>
        private static Modification OGlycan(string residue = "T") =>
            Mod("H1N1", residue, "O-linked glycosylation");

        /// <summary>An O-GlcNAc-like modification: chemically NOT mucin-type, so it must not satisfy the requirement.</summary>
        private static Modification NonGlycan(string residue = "T") => Mod("Phospho", residue, "Test");

        private static Modification Mod(string originalId, string motif, string modificationType)
        {
            Assert.IsTrue(ModificationMotif.TryGetMotif(motif, out ModificationMotif parsedMotif),
                $"test setup: '{motif}' is not a legal modification motif");
            return new Modification(_originalId: originalId, _modificationType: modificationType,
                _target: parsedMotif, _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373);
        }

        /// <summary>
        /// StcE as the enzymology describes it: S/T-X-S/T, cut before the last residue, with an O-glycan
        /// REQUIRED at P2. Registered under a test-only name so the shipped sequence-only StcE entry is
        /// left alone.
        /// </summary>
        private static Protease GlycanAwareStcE(string name)
        {
            var requirement = CleavageRequirement.NonPrime(2, GlycosylationClass.OLinked);
            var motifs = new List<DigestionMotif>
            {
                new("TXT", null, 2, null, requirement),
                new("TXS", null, 2, null, requirement),
                new("SXT", null, 2, null, requirement),
                new("SXS", null, 2, null, requirement),
            };

            var protease = new Protease(name, CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary[name] = protease;
            return protease;
        }

        /// <summary>An OgpA-family protease: cut N-terminal to a glycosylated Ser/Thr, glycan REQUIRED at P1'.</summary>
        private static Protease GlycanAwareOgpA(string name)
        {
            var requirement = CleavageRequirement.Prime(1, GlycosylationClass.OLinked);
            var motifs = new List<DigestionMotif>
            {
                new("T", null, 0, null, requirement),
                new("S", null, 0, null, requirement),
            };

            var protease = new Protease(name, CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary[name] = protease;
            return protease;
        }

        private static List<string> Digest(string sequence, string proteaseName, bool respectPromoting,
            IDictionary<int, List<Modification>> localizedMods = null, List<Modification> variableMods = null,
            InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Retain)
        {
            var protein = localizedMods is null
                ? new Protein(sequence, "TEST")
                : new Protein(sequence, "TEST", oneBasedModifications: localizedMods.ToDictionary(kv => kv.Key, kv => kv.Value));

            var parameters = new DigestionParams(
                protease: proteaseName,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: initiatorMethionineBehavior,
                respectCleavagePromotingModifications: respectPromoting);

            // DISTINCT BASE sequences, not full sequences. A localized modification makes the generator
            // emit the same peptide twice, once carrying it and once not, so counting peptidoforms would
            // conflate "the protease made a second cut" with "the same peptide has two glycoforms" --
            // and those are the two things these tests exist to tell apart.
            return protein.Digest(parameters, new List<Modification>(), variableMods ?? new List<Modification>())
                .Select(p => p.BaseSequence)
                .Distinct()
                .OrderBy(s => s, System.StringComparer.Ordinal)
                .ToList();
        }

        [Test]
        public static void WithTheFlagOn_ATruncationProductBoundaryNeedsNoGlycanToJustifyIt()
        {
            // The signal peptide AAAAK ends after Lys5, so FullDigestion yields the mature chain
            // RPPITQSSL as a "chain start" product. Its N-terminus is a processing site from the
            // database, not a StcE cut, so no StcE motif fits there -- and it must not be asked to.
            // Without a glycan StcE cuts nothing inside it, so the mature chain is the product.
            GlycanAwareStcE("StcE-req-truncation");
            var protein = new Protein("AAAAKRPPITQSSL", "TRUNC",
                proteolysisProducts: new List<TruncationProduct> { new(6, 14, "chain") });
            var parameters = new DigestionParams(protease: "StcE-req-truncation", maxMissedCleavages: 0,
                minPeptideLength: 1, respectCleavagePromotingModifications: true);

            List<string> products = protein.Digest(parameters, new List<Modification>(), new List<Modification>())
                .Select(p => p.BaseSequence).Distinct().ToList();

            CollectionAssert.Contains(products, "RPPITQSSL",
                "the truncation-derived mature chain must survive; got " + string.Join(" | ", products));
        }

        // ---------------------------------------------------------------------------------------
        // The flag is off by default, and off means historically identical
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void TheFlagDefaultsToOff()
        {
            Assert.IsFalse(new DigestionParams().RespectCleavagePromotingModifications,
                "the flag must default to off, so that adding this feature changes no existing search");
        }

        [Test]
        public static void WithTheFlagOff_ARequiringProteaseDigestsExactlyAsTheSequenceMotifAlone()
        {
            GlycanAwareStcE("StcE-req-off");

            // No glycan anywhere, so the real enzyme would not cut at all -- but with the flag off the
            // sequence motif governs, and it does, exactly as the shipped StcE entry does today.
            List<string> products = Digest("RPPITQSSL", "StcE-req-off", respectPromoting: false);

            Assert.AreEqual(2, products.Count, "flag off must reproduce the modification-blind digestion");
        }

        // ---------------------------------------------------------------------------------------
        // The requirement gates the cut: StcE at P2
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void WithTheFlagOn_TheCutIsDroppedWhenTheRequiredGlycanIsAbsent()
        {
            GlycanAwareStcE("StcE-req-absent");

            List<string> products = Digest("RPPITQSSL", "StcE-req-absent", respectPromoting: true);

            // STCE-03 in the truth set: "Deglycosylation abolishes activity."
            Assert.AreEqual(1, products.Count,
                "with no glycan at P2 the protease could not have cut, so the only product is the whole peptide");
            Assert.AreEqual("RPPITQSSL", products[0]);
        }

        [Test]
        public static void WithTheFlagOn_TheCutSurvivesWhenTheRequiredGlycanIsPresentAtP2()
        {
            GlycanAwareStcE("StcE-req-present");

            // RPPITQSSL: the motif is T5-Q6-S7, so P2 is Thr5 and the bond falls before Ser7.
            var mods = new Dictionary<int, List<Modification>> { { 5, new List<Modification> { OGlycan("T") } } };

            List<string> products = Digest("RPPITQSSL", "StcE-req-present", respectPromoting: true, mods);

            // STCE-01: "converting RPPIT*QSSL to RPPIT*Q" -- the two cut products.
            CollectionAssert.Contains(products, "RPPITQ",
                "expected the N-terminal product RPPIT*Q; got " + string.Join(" | ", products));
            CollectionAssert.Contains(products, "SSL",
                "expected the C-terminal product SSL; got " + string.Join(" | ", products));

            // And the READ-THROUGH, which is the third product and the point of slice 4. A localized
            // modification is a site that MAY be occupied, so the digest describes a mixed population:
            // molecules carrying the glycan are cut, molecules without it cannot be, and the intact
            // peptide from the second population is a real product. This assertion counts BASE sequences,
            // so it is not the glycoform of RPPITQ being counted twice -- it is a distinct backbone that
            // spans the site. Before slice 4 it was missing, because at MaxMissedCleavages = 0 the span
            // was never enumerated and the occupancy drop had nothing to put in its place.
            CollectionAssert.Contains(products, "RPPITQSSL",
                "the unglycosylated population cannot be cut, so the intact peptide must survive; got "
                + string.Join(" | ", products));
            Assert.AreEqual(3, products.Count,
                "exactly the two cut products and the read-through; got " + string.Join(" | ", products));
        }

        [Test]
        public static void TheReadThroughCarriesNoMissedCleavage_AndTheGlycosylatedIntactFormDoesNotSurvive()
        {
            // Truth-set STCE-08, and the sharpest statement of what slice 4 does.
            //
            // The read-through is not a missed cleavage. Skipping a site the protease COULD NOT have cut
            // is not a cleavage it missed, so RPPITQSSL comes back at MaxMissedCleavages = 0 with a
            // reported count of zero. That is what lets it exist at all at the caller's budget, and it is
            // why the generation slack that bought its span cannot leak out as an over-budget peptide.
            //
            // The converse matters just as much: the peptidoform that carries the glycan AND spans the
            // site must NOT survive, because with the glycan present StcE would have cut. Keeping it
            // would turn the correction into a no-op that merely added peptides.
            GlycanAwareStcE("StcE-req-readthrough");

            var mods = new Dictionary<int, List<Modification>> { { 5, new List<Modification> { OGlycan("T") } } };
            var protein = new Protein("RPPITQSSL", "TEST",
                oneBasedModifications: mods.ToDictionary(kv => kv.Key, kv => kv.Value));

            var parameters = new DigestionParams(protease: "StcE-req-readthrough", maxMissedCleavages: 0,
                minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                respectCleavagePromotingModifications: true);

            var products = protein.Digest(parameters, new List<Modification>(), new List<Modification>()).ToList();

            var readThrough = products.SingleOrDefault(x => x.BaseSequence == "RPPITQSSL");
            Assert.IsNotNull(readThrough, "the read-through must be generated; got "
                + string.Join(" | ", products.Select(x => x.FullSequence)));
            Assert.AreEqual(0, readThrough.MissedCleavages,
                "an unoccupied site is not a cleavage the protease missed, so the count must be discounted to zero");
            Assert.IsFalse(readThrough.AllModsOneIsNterminus.ContainsKey(5 - 1 + 2),
                "the surviving read-through is the UNGLYCOSYLATED form; the glycosylated one would have been cut");

            Assert.AreEqual(3, products.Select(x => x.BaseSequence).Distinct().Count(),
                "two cut products and one read-through, no more: "
                + string.Join(" | ", products.Select(x => x.FullSequence)));
        }

        [Test]
        public static void WithTheFlagOn_AModificationOfTheWrongClassDoesNotJustifyTheCut()
        {
            GlycanAwareStcE("StcE-req-wrongclass");

            // STCE-02: the same backbone carrying a non-mucin-type modification is NOT cleaved. This is
            // the sharpest case in the truth set -- identical sequence, identical position, one sugar
            // swapped -- and it is the one a sequence-only motif cannot possibly get right.
            var mods = new Dictionary<int, List<Modification>> { { 5, new List<Modification> { NonGlycan("T") } } };

            List<string> products = Digest("RPPITQSSL", "StcE-req-wrongclass", respectPromoting: true, mods);

            Assert.AreEqual(1, products.Count,
                "a modification that is not of the required glycosylation class must not justify the cut");
        }

        [Test]
        public static void WithTheFlagOn_AGlycanAtTheWrongSubsiteDoesNotJustifyTheCut()
        {
            GlycanAwareStcE("StcE-req-wrongsite");

            // Glycan on Ser7, which is P1' of the T5-Q6|S7 cut, not P2. P1' is permitted but never
            // required, and crucially it cannot stand in for the P2 requirement.
            var mods = new Dictionary<int, List<Modification>> { { 7, new List<Modification> { OGlycan("S") } } };

            List<string> products = Digest("RPPITQSSL", "StcE-req-wrongsite", respectPromoting: true, mods);

            Assert.AreEqual(1, products.Count,
                "the requirement names P2; a glycan at P1' is not a substitute, which is exactly the "
                + "distinction a sequence motif cannot express");
        }

        // ---------------------------------------------------------------------------------------
        // The requirement gates the cut: OgpA family at P1'
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void APrimeSideRequirementIsDischargedByThePeptideStartingAtTheCut()
        {
            GlycanAwareOgpA("OgpA-req-present");

            // AHGVTSAPDTRK with the glycan on Thr5: OpeRATOR cuts N-terminal to it, so the products are
            // AHGV and T*SAPDTRK. The constrained residue is the FIRST residue of the second product --
            // the opposite side from StcE, which is the whole reason the requirement carries a side.
            var mods = new Dictionary<int, List<Modification>> { { 5, new List<Modification> { OGlycan("T") } } };

            List<string> products = Digest("AHGVTSAPDTRK", "OgpA-req-present", respectPromoting: true, mods);

            Assert.IsTrue(products.Any(p => p.StartsWith("AHGV", System.StringComparison.Ordinal)),
                "expected the N-terminal product AHGV; got " + string.Join(" | ", products));
            Assert.Greater(products.Count, 1, "the glycosylated Thr justifies the cut before it");
        }

        [Test]
        public static void APrimeSideRequirementDropsEveryCutWhenNoGlycanIsPresent()
        {
            GlycanAwareOgpA("OgpA-req-absent");

            List<string> products = Digest("AHGVTSAPDTRK", "OgpA-req-absent", respectPromoting: true);

            // Every Ser and Thr is a sequence match, so the modification-blind motif shatters this
            // peptide; with the requirement honoured, none of those cuts is justified.
            Assert.AreEqual(1, products.Count,
                "no glycan anywhere means no cut is justified, so the whole peptide survives intact");
            Assert.AreEqual("AHGVTSAPDTRK", products[0]);
        }

        [Test]
        public static void APrimeSideRequirementCutsOnlyWhereTheGlycanIs()
        {
            GlycanAwareOgpA("OgpA-req-selective");

            // The discriminating case for a glycoprotease: the sequence offers four Ser/Thr sites and
            // only one carries a glycan, so exactly one of the four candidate bonds may be severed.
            var mods = new Dictionary<int, List<Modification>> { { 5, new List<Modification> { OGlycan("T") } } };

            List<string> withRequirement = Digest("AHGVTSAPDTRK", "OgpA-req-selective", respectPromoting: true, mods);
            List<string> withoutRequirement = Digest("AHGVTSAPDTRK", "OgpA-req-selective", respectPromoting: false, mods);

            Assert.Less(withRequirement.Count, withoutRequirement.Count,
                "honouring the requirement must produce FEWER peptides than the sequence motif alone -- "
                + "the promoting correction only ever removes peptidoforms, never adds them");
        }

        // ---------------------------------------------------------------------------------------
        // Mixed motifs: a requirement must not leak onto motifs that do not carry one
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void AMotifWithoutARequirementStillCutsWhenAnotherMotifHasOne()
        {
            // StcE-trypsin in miniature: glycan-requiring StcE motifs beside plain tryptic ones. A
            // tryptic cut is answered by the tryptic motif and needs no glycan; only the cuts that
            // nothing but an StcE motif explains have to be justified.
            var requirement = CleavageRequirement.NonPrime(2, GlycosylationClass.OLinked);
            var motifs = new List<DigestionMotif>
            {
                new("TXS", null, 2, null, requirement),
                new("K", null, 1, null),
                new("R", null, 1, null),
            };
            var protease = new Protease("StcE-trypsin-req", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary[protease.Name] = protease;

            // AAKAA: no glycan anywhere, but the tryptic motif matches after Lys3 and carries no
            // requirement, so that cut must still be made.
            List<string> products = Digest("AAKAA", "StcE-trypsin-req", respectPromoting: true);

            Assert.AreEqual(2, products.Count,
                "a motif with no requirement must keep cutting even when a sibling motif has one");
        }

        // ---------------------------------------------------------------------------------------
        // Protein termini are not cuts
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void ProteinTerminiAreNeverTreatedAsUnjustifiedCuts()
        {
            GlycanAwareOgpA("OgpA-req-termini");

            // The full-length peptide's termini are the protein's own, produced by the sequence ending
            // rather than by the protease, so they need no glycan to justify them. If they were checked
            // the whole protein would vanish from the digest.
            List<string> products = Digest("AHGVTSAPDTRK", "OgpA-req-termini", respectPromoting: true);

            Assert.AreEqual(1, products.Count);
            Assert.AreEqual("AHGVTSAPDTRK", products[0],
                "the intact protein must survive: neither of its termini is a protease cut");
        }
        // ---------------------------------------------------------------------------------------
        // Where the glycan may come from, and what happens when it comes from nowhere
        // ---------------------------------------------------------------------------------------

        [Test]
        public static void AVariableModificationMakesASiteFeasibleWithNoDatabaseAnnotation()
        {
            // THE CONFIGURATION THAT BROKE. Feasibility used to be judged only against the database's
            // localized modifications, so a search that supplied the glycan as a VARIABLE modification --
            // the ordinary way to look for O-glycopeptides in MetaMorpheus -- found no annotated glycosite
            // anywhere, judged every site infeasible, filtered the whole site list away, and returned the
            // undigested protein as the only product. The protease was silently switched off.
            //
            // On the parent commit this returns exactly one backbone. Both sources of the modification
            // have to be consulted: the database says "known to be glycosylated here", the configured
            // variable modification says "willing to place a glycan wherever it fits".
            GlycanAwareStcE("StcE-req-varmod");

            List<string> withVariableGlycan = Digest("RPPITQSSL", "StcE-req-varmod", respectPromoting: true,
                variableMods: new List<Modification> { OGlycan() });

            CollectionAssert.Contains(withVariableGlycan, "RPPITQ",
                "a variable O-glycan fits Thr5, so the StcE site is feasible and the cut must be made");
            CollectionAssert.Contains(withVariableGlycan, "SSL",
                "the C-terminal product must survive too: its N-terminal cut is justified by a glycan "
                + "lying in the NEIGHBOURING product, which the out-of-product fallback has to allow for");
        }

        [Test]
        public static void WithNothingAbleToCarryTheGlycanTheGlycoproteaseDoesNotCleave()
        {
            // Pinned deliberately, because the obvious "make it inert" "fix" is wrong and was tried.
            //
            // It is tempting to read the previous test's failure as "the feature must go inert when no
            // glycan is configured" and add a gate mirroring CleavageBlockingPolicy.For
            // AnyConfiguredModificationCanBlockCleavage. That gate breaks the enzymology: StcE, OpeRATOR
            // and IMPa demonstrably do NOT cleave unglycosylated substrate, which is what the truth set's
            // unglycosylated controls (STCE-03, IMPA-12, OGPA-07) encode. The two corrections are not
            // symmetric -- an unconfigured BLOCKING modification cannot remove a site the sequence really
            // has, but an unsatisfiable PROMOTING requirement means there is genuinely no site.
            //
            // So returning almost nothing here is the correct answer, not a defect. What it implies is a
            // PRECONDITION on the flag: a search that turns it on must put the glycan somewhere digestion
            // can see -- annotated in the database or configured as a modification. A glyco search, where
            // the glycan is resolved after identification and never reaches digestion at all, satisfies
            // neither and must not turn the flag on.
            GlycanAwareStcE("StcE-req-noglycan");

            List<string> nothingConfigured = Digest("RPPITQSSL", "StcE-req-noglycan", respectPromoting: true);

            CollectionAssert.AreEqual(new[] { "RPPITQSSL" }, nothingConfigured,
                "with no glycan reachable at all there is no justifiable StcE site, so the substrate "
                + "stays intact -- the published unglycosylated-control result");
        }

        [Test]
        public static void RemovingTheInitiatorMethionineIsNotAProteaseCut()
        {
            // A peptide starting at residue 2 of a Met-initiated sequence got that N-terminus from
            // initiator-methionine removal, not from the protease, so demanding a glycan justify it is a
            // category error. It silently deleted every initiator-cleaved form -- for a glycoprotease,
            // half of all N-terminal peptidoforms -- with no read-through to replace them.
            GlycanAwareStcE("StcE-req-initmet");

            var glycanAtThr6 = new Dictionary<int, List<Modification>> { { 6, new List<Modification> { OGlycan() } } };

            List<string> products = Digest("MRPPITQSSL", "StcE-req-initmet", respectPromoting: true,
                localizedMods: glycanAtThr6, variableMods: new List<Modification> { OGlycan() },
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable);

            CollectionAssert.Contains(products, "RPPITQ",
                "the initiator-cleaved form must survive: its N-terminus is not a protease cut");
            CollectionAssert.Contains(products, "MRPPITQ",
                "and so must the Met-retained form, whose N-terminus is the protein's own");
        }
    }
}
