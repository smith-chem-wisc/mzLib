using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Omics;
using Proteomics;
using Quantification.Differential;

namespace Test.Quantification.Differential;

/// <summary>
/// SETTEST1 S4a (QuantProject <c>DEF-DIFF-SHARE-REGION</c>, <c>-PEPTIDES</c>, <c>-STATE</c>): reading <c>marker_region</c>,
/// the regions' signatures, the marker and reference peptides found by sequence (I and L read as one, ruled by the user
/// 2026-11-14), and a share member's states in the contract's order. Toy databases first, one case per rule; then the
/// share rows of aging's hypotheses-splicing-v2 on an excerpt of the databases their searches read
/// (<c>ReferenceData/PROVENANCE_settest_share.txt</c>), where all 27 named markers must pass.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class SetTestShareTests
{
    // ── toy databases ───────────────────────────────────────────────────────
    //
    // Gene Gx has two forms. C1 holds the exon QQFYSGR at residues 21-27; I1 skips it.
    //   C1: MSTRK AAPEGHK VDLTEWNR QQFYSGR LLNDAMK TTPYWR      (1-5, 6-12, 13-20, 21-27, 28-34, 35-40)
    //   I1: MSTRK AAPEGHK VDLTEWNR LLNDAMK TTPYWR              (I1:20/21 is the exon's skipping junction)
    // O1 (gene Ox) holds TTPYWR and VDITEWNR (VDLTEWNR with I for L); O2 (gene Oy) holds FYSGRLL without the exon's
    // flanks; X has no gene name and holds MSTRK.
    private const string C1 = "MSTRKAAPEGHKVDLTEWNRQQFYSGRLLNDAMKTTPYWR";
    private const string I1 = "MSTRKAAPEGHKVDLTEWNRLLNDAMKTTPYWR";
    private const string O1 = "GGGTTPYWRGGGVDITEWNRGG";
    private const string O2 = "PPPFYSGRLLPPP";
    private const string X = "WWWMSTRKWWW";

    private static ShareSequenceIndex Sequences(bool withLookalike = true) => ShareSequenceIndex.FromSequences(new[]
    {
        ("C1", C1), ("I1", I1), ("O1", withLookalike ? O1 : "GGGTTPYWRGGG"), ("O2", O2), ("X", X),
    });

    private static readonly SetGeneIndex Genes = SetGeneIndex.FromPrimaryNames(new[]
    {
        ("C1", (IEnumerable<string>)new[] { "Gx" }), ("I1", new[] { "Gx" }), ("O1", new[] { "Ox" }), ("O2", new[] { "Oy" }),
        ("X", Array.Empty<string>()),
    });

    private static ShareRegion Range(string accession, int a, int b) =>
        new() { Accession = accession, Kind = ShareRegionKind.Range, First = a, Last = b };

    private static ObservedPeptide P(string sequence) => new(sequence, sequence, new[] { "C1" });

    private static Dictionary<string, SharePeptideClass?> Classes(IEnumerable<string> peptides, string accession, string sequence,
        IReadOnlyList<ShareRegion> regions, ShareSequenceIndex? index = null)
    {
        var observed = peptides.Select(P).ToList();
        var c = SharePeptideClassifier.Classify(observed, sequence, regions, Genes.AccessionsOf("Gx"), "Gx", index ?? Sequences(), Genes);
        return observed.ToDictionary(p => p.FullSequence, p =>
            c.Markers.Contains(p) ? SharePeptideClass.Marker : c.References.Contains(p) ? SharePeptideClass.Reference
            : c.Neither.Contains(p) ? SharePeptideClass.Neither : (SharePeptideClass?)null);
    }

    // ── DEF-DIFF-SHARE-REGION ───────────────────────────────────────────────

    [Test]
    public void EachTokenOfTheAccession_IsARegion_AndEveryOtherWordIsANote()
    {
        var two = ShareRegionReader.Read("Q62261-2:1-36 (alternative N-terminus); Q62261-2:2127-2154 (short C-terminus)", "Q62261-2");
        Assert.That(two.Select(r => (r.Kind, r.First, r.Last)), Is.EqualTo(new[] { (ShareRegionKind.Range, 1, 36), (ShareRegionKind.Range, 2127, 2154) }));

        var junction = ShareRegionReader.Read("P12023-2:junction 288/289 EEVVR|VPTTA (exons 7 and 8 skipped)", "P12023-2").Single();
        Assert.That((junction.Kind, junction.First, junction.Last, junction.Text), Is.EqualTo((ShareRegionKind.Junction, 288, 289, "junction 288/289")));

        var residue = ShareRegionReader.Read("Q8R429:994 (C-terminal residue from exon 22)", "Q8R429").Single();
        Assert.That((residue.Kind, residue.First, residue.Last, residue.Text), Is.EqualTo((ShareRegionKind.Residue, 994, 994, "994")));

        Assert.That(ShareRegionReader.Read("Q8CFN5:87-134 (alpha1 exon; differs from alpha2 Q8CFN5-5)", "Q8CFN5").Single().Text,
            Is.EqualTo("87-134"), "Q8CFN5-5 without a ':' is a note");
        Assert.That(ShareRegionReader.Read("P12023:637-654 (exon 15)", "P12023-2"), Is.Empty, "a token of another accession is a note");
        Assert.That(ShareRegionReader.Read("exon 2 of the long form", "P02688"), Is.Empty);
        Assert.That(ShareRegionReader.Read(null, "P02688"), Is.Empty);
        Assert.That(ShareRegionReader.Read("XP02688:59-84", "P02688"), Is.Empty, "the accession must start the token");
        Assert.That(ShareRegionReader.Read("P02688:99999999999-1", "P02688").Single().LiesIn(1000), Is.False, "a number too large to read");
    }

    [TestCase(ShareRegionKind.Range, 21, 27, 40, true)]
    [TestCase(ShareRegionKind.Range, 27, 21, 40, false)]
    [TestCase(ShareRegionKind.Range, 35, 41, 40, false)]
    [TestCase(ShareRegionKind.Range, 0, 3, 40, false)]
    [TestCase(ShareRegionKind.Residue, 40, 40, 40, true)]
    [TestCase(ShareRegionKind.Junction, 20, 21, 33, true)]
    [TestCase(ShareRegionKind.Junction, 20, 25, 33, false)]
    public void ARegionLiesInItsSequence_OnlyWhereTheContractNumbersIt(ShareRegionKind kind, int a, int b, int length, bool lies)
    {
        var region = new ShareRegion { Accession = "C1", Kind = kind, First = a, Last = b };
        Assert.That(region.LiesIn(length), Is.EqualTo(lies));
    }

    [Test]
    public void ASignature_IsTheRegionWithFiveResiduesEachSide_ClippedAtTheEnds()
    {
        Assert.That(Range("C1", 21, 27).Signature(C1), Is.EqualTo("TEWNRQQFYSGRLLNDA"));
        Assert.That(Range("C1", 2, 3).Signature(C1), Is.EqualTo("MSTRKAAP"), "1 residue before, 5 after");
        Assert.That(new ShareRegion { Accession = "C1", Kind = ShareRegionKind.Residue, First = 40, Last = 40 }.Signature(C1),
            Is.EqualTo("TTPYWR"));
        Assert.That(new ShareRegion { Accession = "I1", Kind = ShareRegionKind.Junction, First = 20, Last = 21 }.Signature(I1),
            Is.EqualTo("TEWNRLLNDA"), "the 5 residues up to a and the 5 from b");
        Assert.That(() => Range("C1", 35, 41).Signature(C1), Throws.ArgumentException.With.Message.Contains("does not lie in"));
    }

    [Test]
    public void APeptideOverlapsARange_ByOneResidue_AndAJunction_ByResidueAOrB()
    {
        var range = Range("C1", 21, 27);
        Assert.That(range.Overlaps(14, 21) && range.Overlaps(27, 30) && !range.Overlaps(13, 20) && !range.Overlaps(28, 34));
        var junction = new ShareRegion { Accession = "I1", Kind = ShareRegionKind.Junction, First = 20, Last = 21 };
        Assert.That(junction.Overlaps(13, 20) && junction.Overlaps(21, 27) && !junction.Overlaps(6, 12) && !junction.Overlaps(22, 27));
    }

    // ── the sequence index ──────────────────────────────────────────────────

    [Test]
    public void ContainmentReadsIAndLAsOneResidue()
    {
        var index = Sequences();
        Assert.That(ShareSequenceIndex.Fold("VDITEWNR"), Is.EqualTo("VDLTEWNR"));
        Assert.That(index.AccessionsContaining("VDLTEWNR"), Is.EqualTo(new[] { "C1", "I1", "O1" }));
        Assert.That(index.AccessionsContaining("VDITEWNR"), Is.EqualTo(new[] { "C1", "I1", "O1" }));
        Assert.That(index.Holds("O1", "VDLTEWNR") && !index.Holds("O2", "VDLTEWNR") && !index.Holds("nope", "A"));
        Assert.That(index.SequenceOf("C1"), Is.EqualTo(C1));
        Assert.That(index.SequenceOf("nope"), Is.Null);
    }

    [Test]
    public void TheIndexFromDatabaseEntries_ReadsTheirBaseSequences()
    {
        var index = ShareSequenceIndex.FromBioPolymers(new IBioPolymer[] { new Protein(C1, "C1"), new Protein(I1, "I1") });
        Assert.That(index.AccessionsContaining("AAPEGHK"), Is.EqualTo(new[] { "C1", "I1" }));
        Assert.That(index.AccessionsContaining("QQFYSGR"), Is.EqualTo(new[] { "C1" }));
    }

    [Test]
    public void AnAccessionTwiceWithOneSequence_IsOneEntry_AndWithTwo_IsRefused()
    {
        var merged = ShareSequenceIndex.FromSequences(new[] { ("P1", "PEPTIDE"), ("P1", "PEPTIDE") });
        Assert.That(merged.AccessionsContaining("PEP"), Is.EqualTo(new[] { "P1" }));
        Assert.That(() => ShareSequenceIndex.FromSequences(new[] { ("P1", "PEPTIDE"), ("P1", "PEPTIDA") }),
            Throws.ArgumentException.With.Message.Contains("two different sequences"));
        Assert.That(() => ShareSequenceIndex.FromSequences(new[] { ("", "PEPTIDE") }), Throws.ArgumentException);
    }

    // ── DEF-DIFF-SHARE-PEPTIDES ─────────────────────────────────────────────

    [Test]
    public void ThePeptidesOfAnExon_AreClassedBySequence()
    {
        var classes = Classes(new[] { "QQFYSGR", "NRQQFY", "FYSGRLL", "AAPEGHK", "LLNDAMK", "VDLTEWNR", "TTPYWR", "MSTRK", "PEGHKVDL", "PPPFY" },
            "C1", C1, new[] { Range("C1", 21, 27) });

        Assert.That(classes["QQFYSGR"], Is.EqualTo(SharePeptideClass.Marker), "only C1 holds it, and C1 holds the exon");
        Assert.That(classes["NRQQFY"], Is.EqualTo(SharePeptideClass.Marker), "overlapping by part is enough");
        Assert.That(classes["FYSGRLL"], Is.EqualTo(SharePeptideClass.Neither), "O2 holds it without the exon's signature");
        Assert.That(classes["AAPEGHK"], Is.EqualTo(SharePeptideClass.Reference));
        Assert.That(classes["LLNDAMK"], Is.EqualTo(SharePeptideClass.Reference));
        Assert.That(classes["PEGHKVDL"], Is.EqualTo(SharePeptideClass.Reference), "no enzyme is assumed");
        Assert.That(classes["VDLTEWNR"], Is.EqualTo(SharePeptideClass.Neither), "O1 holds it with I for L (ruled 2026-11-14)");
        Assert.That(classes["TTPYWR"], Is.EqualTo(SharePeptideClass.Neither), "another gene's entry holds it");
        Assert.That(classes["MSTRK"], Is.EqualTo(SharePeptideClass.Neither), "an entry with no gene name holds it");
        Assert.That(classes["PPPFY"], Is.Null, "not the gene's peptide at all");

        var exact = Classes(new[] { "VDLTEWNR" }, "C1", C1, new[] { Range("C1", 21, 27) }, Sequences(withLookalike: false));
        Assert.That(exact["VDLTEWNR"], Is.EqualTo(SharePeptideClass.Reference), "without the I/L lookalike it is a reference");
    }

    [Test]
    public void APeptideOnlyInSomeOfTheGenesForms_IsNotAReference()
    {
        // Against I1's numbering the exon is gone, so QQFYSGR (C1 only) overlaps nothing and is in only one of the two forms.
        var classes = Classes(new[] { "QQFYSGR", "AAPEGHK" }, "I1", I1, new[] { new ShareRegion { Accession = "I1", Kind = ShareRegionKind.Junction, First = 20, Last = 21 } });
        Assert.That(classes["QQFYSGR"], Is.EqualTo(SharePeptideClass.Neither));
        Assert.That(classes["AAPEGHK"], Is.EqualTo(SharePeptideClass.Reference));
    }

    /// <summary>
    /// H8's case: the skipping junction J2:12/13 joins VDLTEWNR to VPTTAAK, while the long form J1 reads LPTTAAK after its
    /// exon. VPTTAAK starts at residue b and covers only b, yet no other form holds it, so it is a marker; a rule needing both
    /// sides of the junction would lose it. NRVPTT covers a and b.
    /// </summary>
    [Test]
    public void AJunctionMarker_MayCoverResidueBAlone()
    {
        var index = ShareSequenceIndex.FromSequences(new[] { ("J1", "AAAKVDLTEWNRQQFYSGRLPTTAAK"), ("J2", "AAAKVDLTEWNRVPTTAAK") });
        var genes = SetGeneIndex.FromPrimaryNames(new[] { ("J1", (IEnumerable<string>)new[] { "Jg" }), ("J2", new[] { "Jg" }) });
        var junction = new ShareRegion { Accession = "J2", Kind = ShareRegionKind.Junction, First = 12, Last = 13 };
        var peptides = new[] { "VPTTAAK", "NRVPTT", "LPTTAAK", "AAAKVDL" }.Select(P).ToList();

        var c = SharePeptideClassifier.Classify(peptides, "AAAKVDLTEWNRVPTTAAK", new[] { junction }, genes.AccessionsOf("Jg"), "Jg", index, genes);

        Assert.That(c.MarkerBaseSequences, Is.EqualTo(new[] { "NRVPTT", "VPTTAAK" }));
        Assert.That(c.References.Select(p => p.BaseSequence), Is.EqualTo(new[] { "AAAKVDL" }));
        Assert.That(c.Neither.Select(p => p.BaseSequence), Is.EqualTo(new[] { "LPTTAAK" }), "the long form's own, not in J2");
    }

    [Test]
    public void AMultiGeneEntry_IsAnotherGenesEntry_ButMemberAccessionsDecideAlone()
    {
        var index = ShareSequenceIndex.FromSequences(new[] { ("C1", C1), ("I1", I1), ("M", "AAPEGHK") });
        var genes = SetGeneIndex.FromPrimaryNames(new[]
        {
            ("C1", (IEnumerable<string>)new[] { "Gx" }), ("I1", new[] { "Gx" }), ("M", new[] { "Gx", "Gy" }),
        });
        var regions = new[] { Range("C1", 21, 27) };
        var peptides = new[] { P("AAPEGHK") };

        var byName = SharePeptideClassifier.Classify(peptides, C1, regions, genes.AccessionsOf("Gx"), "Gx", index, genes);
        Assert.That(byName.Neither, Has.Count.EqualTo(1), "M also names gene Gy");
        var byAccession = SharePeptideClassifier.Classify(peptides, C1, regions, new[] { "C1", "I1", "M" }, null, index, genes);
        Assert.That(byAccession.References, Has.Count.EqualTo(1));
    }

    [Test]
    public void ARegionOutsideTheSequence_IsRefusedByTheClassifier()
    {
        Assert.That(() => SharePeptideClassifier.Classify(new[] { P("AAPEGHK") }, C1, new[] { Range("C1", 35, 41) }, new[] { "C1" }, "Gx", Sequences(), Genes),
            Throws.ArgumentException.With.Message.Contains("must lie in"));
    }

    // ── DEF-DIFF-SHARE-STATE ────────────────────────────────────────────────

    private const string Header =
        "hypothesis_id\tkind\tprediction\tmember\tspecies\tdataset\tfactor\ttest_level\treference_level\tcontexts\tmrna_log2" +
        "\tindependence\tpost_hoc\tmeasure\tisoform_accession\tmarker_region\tmarker_peptide\tmember_accessions";

    private static SetDatasetContrast Contrast(string prediction = "up", string measure = "exon_peptide", string accession = "C1",
        string region = "C1:21-27 (exon 3)", string listed = "QQFYSGR", string memberAccessions = "")
    {
        var row = string.Join('\t', "H1", "single_event", prediction, "Gx", "mouse", "D1", "age group", "old", "young", "", "",
            "independent", "no", measure, accession, region, listed, memberAccessions);
        return SetTestFreeze.Read(Encoding.UTF8.GetBytes(Header + "\n" + row + "\n"), "t").Parts.Single().Contrasts.Single();
    }

    private static readonly string[] Old = { "o1", "o2", "o3", "o4", "o5" };
    private static readonly string[] Young = { "y1", "y2", "y3", "y4", "y5" };
    private static readonly SetContrastSetup Setup = new() { NumeratorSamples = Old, DenominatorSamples = Young };

    /// <summary>A table over o1-o5 and y1-y5: each peptide's value in the samples named, 20 everywhere it has one.</summary>
    private static ObservationTable Table(params (string Peptide, string[] Samples)[] peptides)
    {
        var samples = Old.Concat(Young).ToList();
        var runs = samples.Select(s => new ObservationRun(s + ".raw", s)).ToList();
        var rows = peptides.Select(p => new PeptideRunValues(p.Peptide, p.Peptide, new[] { "C1" },
            p.Samples.ToDictionary(s => s + ".raw", _ => new RunValue(Math.Pow(2, 20), ObservationState.Quantified)))).ToList();
        return ObservationTable.Build(runs, runs.Select(r => r.FileName), rows, new[] { new ProteinGroupInfo("C1", 0.001, false, false) },
            QuantBasis.MsmsOnly);
    }

    private static readonly string[] Everywhere = Old.Concat(Young).ToArray();

    private static ShareMemberResult State(ObservationTable table, SetDatasetContrast? contrast = null, string? preset = null,
        SetPrediction prediction = SetPrediction.Up)
    {
        var assignment = new SetMemberAssignment { Member = "Gx", ProteinGroup = preset is null ? "C1" : null, State = preset };
        return ShareMemberStates.Apply(contrast ?? Contrast(), prediction, new[] { assignment }, Setup, table, Sequences(), Genes).Single();
    }

    private static ObservationTable Full() => Table(("QQFYSGR", Everywhere), ("AAPEGHK", Everywhere), ("LLNDAMK", Everywhere));

    [Test]
    public void AMemberWithMarkersAndTwoReferencesOnBothSides_IsTested()
    {
        var r = State(Full());
        Assert.That(r.Assignment.State, Is.EqualTo(SetMemberStateFromRow.Tested));
        Assert.That(r.InEstimate);
        Assert.That(r.Peptides!.MarkerBaseSequences, Is.EqualTo(new[] { "QQFYSGR" }));
        Assert.That(r.Peptides.References.Select(p => p.BaseSequence), Is.EqualTo(new[] { "AAPEGHK", "LLNDAMK" }));
        Assert.That((r.MarkerSamplesNumerator, r.MarkerSamplesDenominator, r.ReferenceSamplesNumerator, r.ReferenceSamplesDenominator,
            r.ReferencesWithValue), Is.EqualTo((5, 5, 5, 5, 2)));
        Assert.That(r.ListedMarkerUsed, Is.True);
        Assert.That(r.Regions.Single().Text, Is.EqualTo("21-27"));
    }

    [Test]
    public void AStateFromTheResolver_StandsFirst()
    {
        var r = State(Full(), preset: SetMemberState.NoGroup);
        Assert.That(r.Assignment.State, Is.EqualTo(SetMemberState.NoGroup));
        Assert.That(r.Peptides, Is.Null);
        Assert.That(r.InEstimate, Is.False);
    }

    [Test]
    public void AnIsoformTheSearchDidNotHold_IsNotSearched_ButAnExonAccessionItDidNotHold_IsUnreadable()
    {
        Assert.That(State(Full(), Contrast(measure: "isoform_peptide", accession: "C1-9", region: "C1-9:21-27")).Assignment.State,
            Is.EqualTo(ShareMemberStates.IsoformNotSearched));
        Assert.That(State(Full(), Contrast(accession: "C1-9", region: "C1-9:21-27")).Assignment.State,
            Is.EqualTo(ShareMemberStates.UnreadableRegion));
        Assert.That(State(Full(), Contrast(accession: "", region: "C1:21-27")).Assignment.State,
            Is.EqualTo(ShareMemberStates.UnreadableRegion), "a row with no isoform_accession");
    }

    [TestCase("exon 3 of the long form")]
    [TestCase("C1:35-41 (past the end)")]
    [TestCase("C1:21-27; C1:junction 30/35")]
    public void ARegionThatCannotBeRead_IsUnreadable(string region)
    {
        var r = State(Full(), Contrast(region: region));
        Assert.That(r.Assignment.State, Is.EqualTo(ShareMemberStates.UnreadableRegion));
        Assert.That(r.Peptides, Is.Null);
    }

    [Test]
    public void MarkersSeenOnOneSideOnly_AreSeenOneSide()
    {
        var seen = State(Table(("QQFYSGR", Old), ("AAPEGHK", Everywhere), ("LLNDAMK", Everywhere)));
        Assert.That(seen.Assignment.State, Is.EqualTo(SetMemberStateFromRow.SeenNumeratorOnly));
        Assert.That(seen.InEstimate, Is.False, "never folded into an effect");
        Assert.That(State(Table(("QQFYSGR", Young), ("AAPEGHK", Everywhere), ("LLNDAMK", Everywhere))).Assignment.State,
            Is.EqualTo(SetMemberStateFromRow.SeenDenominatorOnly));
    }

    [Test]
    public void TooFewMarkerValues_OrNoMarker_IsTooFewMarkers()
    {
        Assert.That(State(Table(("QQFYSGR", new[] { "o1", "y1" }), ("AAPEGHK", Everywhere), ("LLNDAMK", Everywhere))).Assignment.State,
            Is.EqualTo(ShareMemberStates.TooFewMarkers));
        Assert.That(State(Table(("AAPEGHK", Everywhere), ("LLNDAMK", Everywhere))).Assignment.State, Is.EqualTo(ShareMemberStates.TooFewMarkers));
    }

    [Test]
    public void OneReferencePeptide_OrReferenceValuesInOneSampleASide_IsTooFewReference()
    {
        Assert.That(State(Table(("QQFYSGR", Everywhere), ("AAPEGHK", Everywhere))).Assignment.State, Is.EqualTo(ShareMemberStates.TooFewReference));
        Assert.That(State(Table(("QQFYSGR", Everywhere), ("AAPEGHK", new[] { "o1", "y1" }), ("LLNDAMK", new[] { "o1", "y1" }))).Assignment.State,
            Is.EqualTo(ShareMemberStates.TooFewReference));
        Assert.That(State(Table(("QQFYSGR", Everywhere), ("AAPEGHK", Everywhere), ("LLNDAMK", Array.Empty<string>()))).ReferencesWithValue,
            Is.EqualTo(1), "a reference with no value in the contrast does not count");
    }

    /// <summary>Markers in 2 of 5 samples a side: enough for a change, too few for a no-change prediction (fewer than half).</summary>
    [Test]
    public void SparseMarkers_AreIncompleteForANoChangePrediction_Only()
    {
        var table = Table(("QQFYSGR", new[] { "o1", "o2", "y1", "y2" }), ("AAPEGHK", Everywhere), ("LLNDAMK", Everywhere));
        var incomplete = State(table, Contrast("unchanged"), prediction: SetPrediction.Unchanged);
        Assert.That(incomplete.Assignment.State, Is.EqualTo(SetMemberStateFromRow.Incomplete));
        Assert.That(incomplete.InEstimate, Is.True, "it stays in the effect and blocks supported (ruled 2026-11-15, as 2026-10-19)");
        Assert.That(State(table).Assignment.State, Is.EqualTo(SetMemberStateFromRow.Tested));
        var sparseReference = Table(("QQFYSGR", Everywhere), ("AAPEGHK", new[] { "o1", "o2", "y1", "y2" }), ("LLNDAMK", new[] { "o1", "o2", "y1", "y2" }));
        Assert.That(State(sparseReference, Contrast("unchanged"), prediction: SetPrediction.Unchanged).Assignment.State,
            Is.EqualTo(SetMemberStateFromRow.Incomplete));
    }

    /// <summary>Each side's markers and reference are checked on their own: sparse on any one side is incomplete.</summary>
    [TestCase("marker", true)]
    [TestCase("marker", false)]
    [TestCase("reference", true)]
    [TestCase("reference", false)]
    public void SparseOnOneSideOnly_IsIncomplete(string kind, bool numerator)
    {
        string[] sparse = numerator ? new[] { "o1", "o2" }.Concat(Young).ToArray() : Old.Concat(new[] { "y1", "y2" }).ToArray();
        var table = kind == "marker"
            ? Table(("QQFYSGR", sparse), ("AAPEGHK", Everywhere), ("LLNDAMK", Everywhere))
            : Table(("QQFYSGR", Everywhere), ("AAPEGHK", sparse), ("LLNDAMK", sparse));
        Assert.That(State(table, Contrast("unchanged"), prediction: SetPrediction.Unchanged).Assignment.State,
            Is.EqualTo(SetMemberStateFromRow.Incomplete));
    }

    /// <summary>A member found by <c>member_accessions</c> has those entries; another gene's entry is any outside them.</summary>
    [Test]
    public void AMemberFoundByItsAccessions_IsClassedOnThoseEntries()
    {
        var r = State(Full(), Contrast(memberAccessions: "C1;I1"));
        Assert.That(r.Assignment.State, Is.EqualTo(SetMemberStateFromRow.Tested));
        Assert.That(r.Peptides!.References.Select(p => p.BaseSequence), Is.EqualTo(new[] { "AAPEGHK", "LLNDAMK" }));
        Assert.That(State(Full(), Contrast(memberAccessions: "C1")).Peptides!.References.Select(p => p.BaseSequence),
            Is.Empty, "I1, outside C1, also holds them, so they are in another gene's entry");
    }

    [Test]
    public void TheListedMarkers_AreReportedAgainstThoseFound_NeverUsed()
    {
        Assert.That(State(Full(), Contrast(listed: "NOTFOUND")).ListedMarkerUsed, Is.False);
        Assert.That(State(Full(), Contrast(listed: "QQFYSGR;NOTFOUND")).ListedMarkerUsed, Is.False, "every listed marker, not any");
        Assert.That(State(Full(), Contrast(listed: "NOTFOUND")).Assignment.State, Is.EqualTo(SetMemberStateFromRow.Tested));
        Assert.That(State(Full(), Contrast(listed: "")).ListedMarkerUsed, Is.Null);
        Assert.That(Contrast(listed: "QQFYSGR; NRQQFY").Rows.Single().MarkerPeptides, Is.EqualTo(new[] { "QQFYSGR", "NRQQFY" }));
    }

    [Test]
    public void AMemberNotInTheContrast_OrASampleNotInTheTable_IsRefused()
    {
        var stranger = new SetMemberAssignment { Member = "Zz", ProteinGroup = "C1" };
        Assert.That(() => ShareMemberStates.Apply(Contrast(), SetPrediction.Up, new[] { stranger }, Setup, Full(), Sequences(), Genes),
            Throws.ArgumentException.With.Message.Contains("is not a member"));
        var setup = new SetContrastSetup { NumeratorSamples = new[] { "zz" }, DenominatorSamples = Young };
        Assert.That(() => ShareMemberStates.Apply(Contrast(), SetPrediction.Up, Array.Empty<SetMemberAssignment>(), setup, Full(), Sequences(), Genes),
            Throws.ArgumentException.With.Message.Contains("no sample 'zz'"));
    }

    // ── the real rows on the real databases ─────────────────────────────────

    private static string Data(string file) =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "Quantification", "Differential", "ReferenceData", file);

    /// <summary>
    /// hypotheses-splicing-v2's 21 share rows, each on its species' databases as the searches read them: every named marker
    /// overlaps its region and passes the marker rule (27 of 27, v3.15 and v3.16). H8's VPTTAASTPDAVDK starts at residue
    /// 289 of P12023-2 and P08592-2; App exon 15's GLTTRPGSGLTNIK is also in P12023-2 and P12023-3, which hold the exon.
    /// </summary>
    [Test]
    public void EveryNamedMarkerOfTheSplicingFreeze_PassesOnTheSearchedDatabases()
    {
        var freeze = SetTestFreeze.ReadFile(Data("settest_share_rows.tsv"), "hypotheses-splicing-v2");
        var lines = File.ReadAllLines(Data("settest_share_entries.tsv")).Skip(1).Where(l => l.Length > 0).Select(l => l.Split('\t')).ToList();
        var bySpecies = lines.GroupBy(c => c[0]).ToDictionary(g => g.Key, g => (
            Sequences: ShareSequenceIndex.FromSequences(g.Select(c => (c[2], c[4]))),
            Genes: SetGeneIndex.FromPrimaryNames(g.Select(c => (c[2], (IEnumerable<string>)c[3].Split(';', StringSplitOptions.RemoveEmptyEntries))))));

        int named = 0;
        foreach (var row in freeze.Parts.SelectMany(p => p.Contrasts).SelectMany(c => c.Rows))
        {
            var (sequences, genes) = bySpecies[row.Species.Trim().ToLowerInvariant()];
            string accession = row.IsoformAccession!;
            string sequence = sequences.SequenceOf(accession) ?? throw new AssertionException($"{accession} is not in the excerpt");
            var regions = ShareRegionReader.Read(row.MarkerRegion, accession);
            Assert.That(regions, Is.Not.Empty, $"{row.HypothesisId} {row.Member}");
            Assert.That(regions.All(r => r.LiesIn(sequence.Length)), $"{row.HypothesisId} {row.Member} {accession}");
            var peptides = row.MarkerPeptides.Select(m => new ObservedPeptide(m, m, new[] { accession })).ToList();
            var classes = SharePeptideClassifier.Classify(peptides, sequence, regions, genes.AccessionsOf(row.Member), row.Member, sequences, genes);
            Assert.That(classes.MarkerBaseSequences, Is.EqualTo(row.MarkerPeptides.Order(StringComparer.Ordinal)),
                $"{row.HypothesisId} {row.Member} {row.Dataset} {accession}");
            named += row.MarkerPeptides.Count;
        }
        Assert.That(named, Is.EqualTo(27));
        var mouse = bySpecies["mouse"].Sequences;
        Assert.That(ShareRegionReader.Read("P12023-2:junction 288/289 EEVVR|VPTTA", "P12023-2").Single().Signature(mouse.SequenceOf("P12023-2")!),
            Is.EqualTo("EEVVRVPTTA"), "the signature aging writes in the cell");
        Assert.That(freeze.Parts.SelectMany(p => p.Contrasts).SelectMany(c => c.Rows).Count(), Is.EqualTo(21));
    }
}
