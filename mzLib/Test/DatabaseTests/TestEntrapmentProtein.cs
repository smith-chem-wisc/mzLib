using System;
using System.Globalization;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using UsefulProteomicsDatabases.EntrapmentGeneration;

namespace Test.DatabaseTests;

[TestFixture]
[ExcludeFromCodeCoverage]
public class EntrapmentProteinTests
{
    private static readonly HashSet<string> NothingForbidden = new();
    private const string Sequence = "SYKALADQMNLLLSKGGVDTTPFAWENDRQISTLGGYK";

    private static IDigestionParams Tryptic => new DigestionParams("trypsin", minPeptideLength: 7, maxMissedCleavages: 2);

    private static Modification Phospho()
    {
        ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
        return new Modification(_originalId: "Phospho", _modificationType: "Common Biological",
            _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);
    }

    // ---- accession ---------------------------------------------------------

    [Test]
    public void Accession_PutsTheMarkerFirstSoTheLoaderFindsIt()
    {
        string accession = EntrapmentAccession.Format("P12345", fold: 0);

        // Prefix, because ProteinDbLoader detects entrapment by looking for the identifier anywhere
        // in the accession and prepends it when absent -- matching the convention means a written
        // database reloads as entrapment with no extra argument, and is never double-prefixed.
        Assert.That(accession, Does.StartWith(ProteinDbLoader.DefaultEntrapmentIdentifier));
        Assert.That(accession, Does.Contain("P12345"));
    }

    [Test]
    public void Accession_RoundTripsTheTargetAndFold()
    {
        for (int fold = 0; fold < 4; fold++)
        {
            string accession = EntrapmentAccession.Format("P12345", fold);
            Assert.That(EntrapmentAccession.TryParse(accession, out string target, out int parsedFold), Is.True);
            Assert.That(target, Is.EqualTo("P12345"));
            Assert.That(parsedFold, Is.EqualTo(fold));
        }
    }

    [Test]
    public void Accession_RefusesSomethingThatIsNotAnEntrapmentAccession()
    {
        Assert.That(EntrapmentAccession.TryParse("P12345", out _, out _), Is.False);
        Assert.That(EntrapmentAccession.TryParse("", out _, out _), Is.False);
        Assert.That(EntrapmentAccession.TryParse(null, out _, out _), Is.False);
    }

    [Test]
    public void Accession_SurvivesATargetAccessionContainingUnderscores()
    {
        // UniProt entry names always carry one, so splitting naively would lose the fold.
        string accession = EntrapmentAccession.Format("SP_HUMAN_P12345", fold: 3);
        Assert.That(EntrapmentAccession.TryParse(accession, out string target, out int fold), Is.True);
        Assert.That(target, Is.EqualTo("SP_HUMAN_P12345"));
        Assert.That(fold, Is.EqualTo(3));
    }

    [Test]
    public void Accession_RefusesToFormatSomethingUnusable()
    {
        var noTarget = Assert.Throws<MzLibUtil.MzLibException>(() => EntrapmentAccession.Format("", 0));
        Assert.That(noTarget!.Message, Does.Contain("target accession"));

        var negative = Assert.Throws<MzLibUtil.MzLibException>(() => EntrapmentAccession.Format("P12345", -1));
        Assert.That(negative!.Message, Does.Contain("negative"));
    }

    [TestCase("Random_P12345", TestName = "Accession parse - no fold marker")]
    [TestCase("Random_P12345_f", TestName = "Accession parse - fold marker with no digits")]
    [TestCase("Random_P12345_fX", TestName = "Accession parse - fold that is not a number")]
    [TestCase("Random_P12345_f-1", TestName = "Accession parse - negative fold")]
    [TestCase("Random__f1", TestName = "Accession parse - empty target accession")]
    [TestCase("Random_f1", TestName = "Accession parse - marker but nothing before the fold")]
    public void Accession_RefusesMalformedInput(string accession)
    {
        Assert.That(EntrapmentAccession.TryParse(accession, out _, out _), Is.False);
    }

    // ---- protein generation ------------------------------------------------

    [Test]
    public void Create_MakesAnEntrapmentProteinIsomericWithItsTarget()
    {
        var target = new Protein(Sequence, "P12345");
        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.IsEntrapment, Is.True);
        Assert.That(entrapment.BaseSequence, Is.Not.EqualTo(target.BaseSequence));
        Assert.That(string.Concat(entrapment.BaseSequence.OrderBy(c => c)),
            Is.EqualTo(string.Concat(target.BaseSequence.OrderBy(c => c))));
        Assert.That(EntrapmentAccession.TryParse(entrapment.Accession, out string parsed, out _), Is.True);
        Assert.That(parsed, Is.EqualTo("P12345"));
    }

    [Test]
    public void Create_CarriesModificationsToTheResidueTheyWereOn()
    {
        Modification phospho = Phospho();
        var mods = new Dictionary<int, List<Modification>>();
        for (int oneBased = 1; oneBased <= Sequence.Length; oneBased++)
        {
            if (Sequence[oneBased - 1] == 'T')
            {
                mods[oneBased] = new List<Modification> { phospho };
            }
        }
        Assert.That(mods, Is.Not.Empty, "the fixture must actually carry modifications");

        var target = new Protein(Sequence, "P12345", oneBasedModifications: mods);
        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(mods.Count),
            "every modification on a retained residue must survive the move");
        foreach (int position in entrapment.OneBasedPossibleLocalizedModifications.Keys)
        {
            Assert.That(entrapment.BaseSequence[position - 1], Is.EqualTo('T'),
                "a modification must land on the residue its motif requires");
        }
    }

    [Test]
    public void Create_DropsModificationsOnExcisedResidues()
    {
        // "SSSSSSR" has no partner and is long enough to be identified, so it is excised -- and any
        // modification sitting on it goes with it rather than landing somewhere arbitrary.
        const string withHomopolymer = "SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR";
        ModificationMotif.TryGetMotif("S", out ModificationMotif motif);
        var onExcisedTract = new Modification(_originalId: "Phospho", _modificationType: "Common Biological",
            _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 79.96633);

        var target = new Protein(withHomopolymer, "P12345", oneBasedModifications:
            new Dictionary<int, List<Modification>> { { 17, new List<Modification> { onExcisedTract } } });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.BaseSequence, Does.Not.Contain("SSSSSSR"));
        Assert.That(entrapment.OneBasedPossibleLocalizedModifications, Is.Empty);
    }

    [Test]
    public void Create_GivesEachFoldItsOwnAccessionAndSequence()
    {
        var target = new Protein(Sequence, "P12345");

        var entrapments = Enumerable.Range(0, 4)
            .Select(k => EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden, fold: k, foldCount: 4))
            .ToList();

        Assert.That(entrapments.Select(p => p.Accession).Distinct().Count(), Is.EqualTo(4));
        Assert.That(entrapments.Select(p => p.BaseSequence).Distinct().Count(), Is.EqualTo(4));
    }

    [Test]
    public void Create_RefusesANullTarget()
    {
        var thrown = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentProteinGenerator.Create(null, Tryptic, NothingForbidden));
        Assert.That(thrown!.Message, Does.Contain("null target"));
    }

    [Test]
    public void Create_HandlesAProteinCarryingNoModifications()
    {
        var target = new Protein(Sequence, "P12345");
        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);
        Assert.That(entrapment.OneBasedPossibleLocalizedModifications, Is.Empty);
    }

    [Test]
    public void Create_KeepsAModificationThatLandsOnTheVeryFirstResidue()
    {
        // Position zero is a real destination, not a "missing" marker -- only a negative means
        // excised. Conflating them would silently drop a modification whenever its residue moved to
        // the front. "AK" is below the minimum length and has no rearrangement, so it is kept
        // verbatim with an identity map: its first residue is GUARANTEED to land at position zero,
        // which an ordinary rearranged piece would only reach by luck.
        ModificationMotif.TryGetMotif("A", out ModificationMotif motif);
        var onFirstResidue = new Modification(_originalId: "Acetyl", _modificationType: "Common Biological",
            _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 42.010565);

        var target = new Protein("AKGGVDTTPFAWENDRQISTLGGYK", "P12345", oneBasedModifications:
            new Dictionary<int, List<Modification>> { { 1, new List<Modification> { onFirstResidue } } });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden,
            out EntrapmentAssembly assembly);

        Assert.That(assembly.TargetToEntrapmentPosition[0], Is.Zero, "the fixture must exercise position zero");
        Assert.That(entrapment.OneBasedPossibleLocalizedModifications.ContainsKey(1), Is.True,
            "a modification landing at position zero must survive");
        Assert.That(entrapment.BaseSequence[0], Is.EqualTo('A'));
    }

    [Test]
    public void Create_DropsPositionalAnnotationsThatTheRearrangementInvalidates()
    {
        // Sequence variations, truncation products, disulfide bonds and splice sites all describe
        // the TARGET's sequence. Once the residues move they mean nothing, and once excision
        // shortens the protein a coordinate can point PAST ITS END -- MetaMorpheus applies sequence
        // variations while loading a database, indexes with that coordinate and throws. Measured on
        // the human proteome: 131 entrapment entries carried a feature position beyond their own
        // length, and every search against that database failed to load it.
        var target = new Protein("SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR", "P12345",
            proteolysisProducts: new List<TruncationProduct> { new(1, 36, "full-length") },
            disulfideBonds: new List<DisulfideBond> { new(3, 30, "bond") },
            spliceSites: new List<SpliceSite> { new(5, 9, "site") });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.TruncationProducts, Is.Empty);
        Assert.That(entrapment.DisulfideBonds, Is.Empty);
        Assert.That(entrapment.SpliceSites, Is.Empty);
        Assert.That(entrapment.SequenceVariations, Is.Empty);
        Assert.That(entrapment.AppliedSequenceVariations, Is.Empty);

        // The target keeps everything -- we copy, we do not mutate.
        Assert.That(target.TruncationProducts, Is.Not.Empty);
        Assert.That(target.DisulfideBonds, Is.Not.Empty);
    }

    [Test]
    public void Create_NeverLeavesAnAnnotationPointingPastTheSequence()
    {
        // The failure mode directly: the fixture excises SSSSSSR, so the entrapment protein is
        // seven residues shorter than the annotations that described its target.
        var target = new Protein("SYKALADQMNLLLSKSSSSSSRGGVDTTPFAWENDR", "P12345",
            proteolysisProducts: new List<TruncationProduct> { new(30, 36, "C-terminal chunk") });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.BaseSequence.Length, Is.LessThan(target.BaseSequence.Length));

        // Assert the contract directly rather than looping a collection the fix empties. The old
        // form iterated TruncationProducts asserting each end position was in range, and this path
        // clears that collection -- so the loop body never ran and the only live assertion was that
        // the sequence got shorter. A test named for the fix could not observe it.
        Assert.That(entrapment.TruncationProducts, Is.Empty,
            "positional annotations describe the target's coordinates and are dropped, not remapped");
        Assert.That(entrapment.SequenceVariations, Is.Empty);
        Assert.That(entrapment.DisulfideBonds, Is.Empty);
        Assert.That(entrapment.SpliceSites, Is.Empty);

        // And the property that mattered, stated over everything the entry actually carries.
        IEnumerable<int> everyPosition = entrapment.TruncationProducts
            .Where(t => t.OneBasedEndPosition.HasValue).Select(t => t.OneBasedEndPosition.Value)
            .Concat(entrapment.DisulfideBonds.Select(b => b.OneBasedEndPosition))
            .Concat(entrapment.SpliceSites.Select(x => x.OneBasedEndPosition))
            .Concat(entrapment.OneBasedPossibleLocalizedModifications.Keys);
        Assert.That(everyPosition, Is.All.LessThanOrEqualTo(entrapment.BaseSequence.Length));
    }

    [Test]
    public void Proteoform_CarriesTruncationProductsAndKeepsThemInRange()
    {
        // The same property where the collection is NOT empty, so the range assertion has something
        // to bite on. Proteoform mode excises nothing, so the length never changes and the spans
        // stay valid -- which is why that path can carry them at all.
        var target = new Protein("MADCQVLGYTTPDNRAWEDSFLCGKQPTMLNDVAHERGGLYT", "P12345",
            proteolysisProducts: new List<TruncationProduct> { new(2, 25, "chain") });

        Protein entrapment = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden, out _);

        Assert.That(entrapment.TruncationProducts, Is.Not.Empty, "or the assertion below is vacuous");
        Assert.That(entrapment.TruncationProducts.Select(t => t.OneBasedEndPosition!.Value),
            Is.All.LessThanOrEqualTo(entrapment.BaseSequence.Length));
    }

    // ---- terminal modifications --------------------------------------------

    private static Modification TerminalMod(string motifResidue, string restriction, string id) 
    {
        ModificationMotif.TryGetMotif(motifResidue, out ModificationMotif motif);
        return new Modification(_originalId: id, _modificationType: "Common Biological",
            _target: motif, _locationRestriction: restriction, _monoisotopicMass: 42.010565);
    }

    [Test]
    public void Create_KeepsAModificationRestrictedToTheProteinNTerminus()
    {
        // A rearrangement that moved this residue off the terminus would make the modification
        // invalid for its location, and mzLib would drop it -- silently. Measured before anchoring:
        // 3,946 N-terminal modifications lost across the human proteome, 2.4% of all of them.
        var target = new Protein("MAAALGGDRKGGVDTTPFAWENDRQISTLGGYK", "P12345",
            oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { TerminalMod("M", "N-terminal.", "N-acetylmethionine") } }
            });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.OneBasedPossibleLocalizedModifications.ContainsKey(1), Is.True,
            "an N-terminally restricted modification must still be at position 1");
        Assert.That(entrapment.BaseSequence[0], Is.EqualTo('M'));
    }

    [Test]
    public void Create_KeepsAModificationOnTheSecondResidue()
    {
        // Position 2 matters as much as position 1: a modification annotated after initiator
        // methionine cleavage lands on the second residue. Both are anchored.
        //
        // The second residue must not be one of a run of identical residues, or the test passes
        // whether or not anchoring works -- with "MAAA..." the unranking returns an A to slot 1
        // anyway and the assertion cannot fail. Here the S is the only one in its piece, so it stays
        // at position 2 only because it is held there.
        var target = new Protein("MSAALGGDRKGGVDTTPFAWENDRQIATLGGYK", "P12345",
            oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 2, new List<Modification> { TerminalMod("S", "N-terminal.", "N-acetylserine") } }
            });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.OneBasedPossibleLocalizedModifications.ContainsKey(2), Is.True);
        Assert.That(entrapment.BaseSequence[1], Is.EqualTo('S'));
        Assert.That(entrapment.BaseSequence.Count(c => c == 'S'), Is.EqualTo(1),
            "the fixture only discriminates while its second residue is unique in the sequence");
    }

    [Test]
    public void Create_KeepsAModificationRestrictedToTheProteinCTerminus()
    {
        const string sequence = "MAAALGGDRKGGVDTTPFAWENDRQISTLGGYA";
        var target = new Protein(sequence, "P12345",
            oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { sequence.Length, new List<Modification> { TerminalMod("A", "C-terminal.", "Amidation") } }
            });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.OneBasedPossibleLocalizedModifications.ContainsKey(entrapment.BaseSequence.Length),
            Is.True, "a C-terminally restricted modification must still be at the last position");
        Assert.That(entrapment.BaseSequence[^1], Is.EqualTo('A'));
    }

    [Test]
    public void Create_ConservesTheModificationCountWhenTerminalModificationsArePresent()
    {
        // The assertion the plan asks for by name: not "the shortfall is reported", but conserved.
        const string sequence = "MAAALGGDRKGGVDTTPFAWENDRQISTLGGYK";
        var target = new Protein(sequence, "P12345",
            oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { TerminalMod("M", "N-terminal.", "N-acetylmethionine") } },
                { 2, new List<Modification> { TerminalMod("A", "N-terminal.", "N-acetylalanine") } },
                { 18, new List<Modification> { Phospho() } }
            });
        int before = target.OneBasedPossibleLocalizedModifications.Sum(kv => kv.Value.Count);

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(entrapment.OneBasedPossibleLocalizedModifications.Sum(kv => kv.Value.Count),
            Is.EqualTo(before), "no modification may be lost when none of them sits on an excised residue");
    }

    [Test]
    public void Create_AnchorsTheTerminiWhetherOrNotAnythingIsModified()
    {
        // Anchoring is unconditional. Pinning only where a terminal modification happens to sit
        // would make the entrapment sequence a function of the annotations as well as the residues,
        // so two databases over the same proteome with different annotations would disagree on
        // their sequences -- and the determinism the pairing rests on would quietly weaken.
        const string sequence = "MAAALGGDRKGGVDTTPFAWENDRQISTLGGYK";
        Protein bare = EntrapmentProteinGenerator.Create(new Protein(sequence, "P1"), Tryptic, NothingForbidden);
        Protein annotated = EntrapmentProteinGenerator.Create(
            new Protein(sequence, "P1", oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { TerminalMod("M", "N-terminal.", "N-acetylmethionine") } }
            }), Tryptic, NothingForbidden);

        Assert.That(annotated.BaseSequence, Is.EqualTo(bare.BaseSequence),
            "the sequence must not depend on whether the protein carried a terminal modification");
        Assert.That(bare.BaseSequence[0], Is.EqualTo('M'));
        Assert.That(bare.BaseSequence[1], Is.EqualTo('A'));
        Assert.That(bare.BaseSequence[^1], Is.EqualTo(sequence[^1]));
    }

    [Test]
    public void Create_AnchorsBothEndsOfASingleBasePieceProtein()
    {
        // A protein with no internal cleavage site is its own first AND last piece, so both
        // anchoring branches fire at once. Nothing else in the fixtures exercises that.
        const string sequence = "MAAALGGDQISTLGGYA";   // no K or R -- one piece
        var target = new Protein(sequence, "P12345",
            oneBasedModifications: new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { TerminalMod("M", "N-terminal.", "N-acetylmethionine") } },
                { sequence.Length, new List<Modification> { TerminalMod("A", "C-terminal.", "Amidation") } }
            });

        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden,
            out EntrapmentAssembly assembly);

        Assert.That(assembly.Pieces.Count, Is.EqualTo(1), "the fixture must actually be a single piece");
        Assert.That(entrapment.BaseSequence[0], Is.EqualTo('M'));
        Assert.That(entrapment.BaseSequence[1], Is.EqualTo('A'));
        Assert.That(entrapment.BaseSequence[^1], Is.EqualTo('A'));
        Assert.That(entrapment.BaseSequence, Is.Not.EqualTo(sequence), "the interior must still move");
        Assert.That(entrapment.OneBasedPossibleLocalizedModifications.Sum(kv => kv.Value.Count),
            Is.EqualTo(2), "both terminal modifications must survive");
    }

    [Test]
    public void Create_AnchoringActuallyShrinksThePermutationSpace()
    {
        // If anchoring silently did nothing, every test above could still pass by luck on a short
        // fixture. Hold the termini and the space must be strictly smaller.
        var motifs = global::Omics.Digestion.DigestionMotif.ParseDigestionMotifsFromString("K|,R|");
        const string piece = "MAAALGGDQISTLGGYA";

        var free = global::UsefulProteomicsDatabases.DecoySequenceValidator.PermutationSpaceSize(piece, motifs);
        var held = global::UsefulProteomicsDatabases.DecoySequenceValidator.PermutationSpaceSize(piece, motifs,
            new[] { 0, 1, piece.Length - 1 });

        Assert.That(held, Is.LessThan(free));
        Assert.That(held, Is.GreaterThan(System.Numerics.BigInteger.One));
    }

    // ---- pairing -----------------------------------------------------------

    [Test]
    public void Pairing_RecoversTheTargetPeptideFromCompositionAlone()
    {
        var target = new Protein(Sequence, "P12345");
        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        var pairing = new EntrapmentPairing(target, Tryptic);

        foreach (var peptide in entrapment.Digest(Tryptic, new List<Modification>(), new List<Modification>()))
        {
            string entrapmentPeptide = peptide.BaseSequence;
            Assert.That(pairing.TryResolve(entrapmentPeptide, out string targetPeptide), Is.True,
                "'" + entrapmentPeptide + "' should resolve to the target peptide it was built from");
            Assert.That(string.Concat(targetPeptide.OrderBy(c => c)),
                Is.EqualTo(string.Concat(entrapmentPeptide.OrderBy(c => c))));
        }
    }

    [Test]
    public void Pairing_NeedsNoSideFileBeyondTheTargetDatabase()
    {
        // The whole point of pairing on composition: a search reports a protein accession and a
        // peptide sequence, and that is enough. Nothing has to be carried alongside the database.
        var target = new Protein(Sequence, "P12345");
        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        Assert.That(EntrapmentAccession.TryParse(entrapment.Accession, out string targetAccession, out _), Is.True);
        Assert.That(targetAccession, Is.EqualTo(target.Accession));

        var pairing = new EntrapmentPairing(target, Tryptic);
        string first = entrapment.Digest(Tryptic, new List<Modification>(), new List<Modification>())
            .First().BaseSequence;
        Assert.That(pairing.TryResolve(first, out _), Is.True);
    }

    [Test]
    public void Pairing_FlagsPeptidesOfIdenticalCompositionAsAmbiguous()
    {
        // P57071 really does contain LIHTGVK and LIHTVGK: same residues, same pinned K, different
        // order. Composition cannot tell them apart, so neither is paired -- they are reported.
        var target = new Protein("LIHTGVKAAADEFGHIKLIHTVGKMMWWYYPPQQR", "P57071");
        var pairing = new EntrapmentPairing(target, Tryptic);

        Assert.That(pairing.AmbiguousPeptides, Does.Contain("LIHTGVK"));
        Assert.That(pairing.AmbiguousPeptides, Does.Contain("LIHTVGK"));
        Assert.That(pairing.TryResolve("LIHTVGK", out _), Is.False,
            "an ambiguous key must refuse to guess");
    }

    [Test]
    public void Pairing_RefusesANullTarget()
    {
        var thrown = Assert.Throws<MzLibUtil.MzLibException>(() => new EntrapmentPairing(null, Tryptic));
        Assert.That(thrown!.Message, Does.Contain("target protein"));
    }

    [Test]
    public void Pairing_RefusesAPeptideItHasNeverSeen()
    {
        var pairing = new EntrapmentPairing(new Protein(Sequence, "P12345"), Tryptic);
        Assert.That(pairing.TryResolve("WWWWWWWW", out _), Is.False);
        Assert.That(pairing.TryResolve("", out _), Is.False);
    }

    // ---- round trip --------------------------------------------------------

    [Test]
    public void EntrapmentProtein_ReloadsAsEntrapmentFromXml()
    {
        var target = new Protein(Sequence, "P12345");
        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "entrapment_roundtrip.xml");
        ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
            new List<Protein> { target, entrapment }, path);

        List<Protein> reloaded = ProteinDbLoader.LoadProteinXML(path, true, DecoyType.None,
            new List<Modification>(), false, new List<string>(), out _);

        Protein reloadedEntrapment = reloaded.Single(p => p.Accession != "P12345");
        Assert.That(reloadedEntrapment.IsEntrapment, Is.True,
            "the accession alone must be enough for the loader to know this is entrapment");
        Assert.That(reloadedEntrapment.BaseSequence, Is.EqualTo(entrapment.BaseSequence));
        Assert.That(EntrapmentAccession.TryParse(reloadedEntrapment.Accession, out string parsed, out _), Is.True);
        Assert.That(parsed, Is.EqualTo("P12345"));

        File.Delete(path);
    }

    [Test]
    public void EntrapmentProtein_PairingSurvivesTheRoundTrip()
    {
        var target = new Protein(Sequence, "P12345");
        Protein entrapment = EntrapmentProteinGenerator.Create(target, Tryptic, NothingForbidden);

        string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "entrapment_pairing.xml");
        ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
            new List<Protein> { target, entrapment }, path);
        List<Protein> reloaded = ProteinDbLoader.LoadProteinXML(path, true, DecoyType.None,
            new List<Modification>(), false, new List<string>(), out _);

        Protein reloadedTarget = reloaded.Single(p => p.Accession == "P12345");
        Protein reloadedEntrapment = reloaded.Single(p => p.Accession != "P12345");

        // Write, read, digest, pair -- with nothing carried alongside the file.
        var pairing = new EntrapmentPairing(reloadedTarget, Tryptic);
        foreach (var peptide in reloadedEntrapment.Digest(Tryptic, new List<Modification>(), new List<Modification>()))
        {
            Assert.That(pairing.TryResolve(peptide.BaseSequence, out string targetPeptide), Is.True);
            Assert.That(string.Concat(targetPeptide.OrderBy(c => c)),
                Is.EqualTo(string.Concat(peptide.BaseSequence.OrderBy(c => c))));
        }

        File.Delete(path);
    }

    /// <summary>
    /// A protein carrying sequence variants, expanded the way LoadProteinXML expands one: several
    /// proteins sharing a single consensus. This is the shape a real XML database has, and the
    /// shape a FASTA database never has -- which is why the defect below only ever showed on the
    /// XML arm.
    /// </summary>
    private static List<Protein> ExpandedVariantProteins()
    {
        var insertion = new SequenceVariation(5, 5, "A", "AAA", "Insertion at 5");
        var substitution = new SequenceVariation(10, 10, "G", "R", "Substitution at 10");
        var consensus = new Protein("MAAAAGAAAAGKPEPTIDESAMPLERTESTPEPTIDEK", "P00001",
            sequenceVariations: new List<SequenceVariation> { insertion, substitution });

        return VariantApplication
            .ApplyAllVariantCombinations(consensus, new List<SequenceVariation> { insertion, substitution },
                maxCombinations: 10)
            .Cast<Protein>()
            .ToList();
    }

    [Test]
    public void ExpandedVariantsShareOneConsensus()
    {
        // Guard on the fixture itself. A fixture that quietly stopped expanding would make every
        // assertion below pass without exercising anything -- which is how a previous test in this
        // suite came to cover a case it never reached.
        List<Protein> targets = ExpandedVariantProteins();

        Assert.That(targets.Count, Is.GreaterThan(1),
            "fixture must expand into several proteins, or it does not exercise the case");
        Assert.That(targets.Select(p => p.ConsensusVariant).Distinct().Count(), Is.EqualTo(1),
            "the expanded proteins must share a single consensus");
    }

    [Test]
    public void OneEntrapmentProteinPerEntryNotPerAppliedVariant()
    {
        List<Protein> targets = ExpandedVariantProteins();

        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateEntrapment(
            targets, Tryptic, NothingForbidden);

        // One per database ENTRY. The target list holds several proteins but the database holds one
        // entry, because variants are annotations on a consensus rather than entries of their own.
        Assert.That(entrapment.Count, Is.EqualTo(1));
        Assert.That(entrapment[0].IsEntrapment, Is.True);
        Assert.That(entrapment[0].Accession, Is.EqualTo("Random_P00001_f0"),
            "the accession must name the consensus, not an applied variant");
    }

    [Test]
    public void EntrapmentEntriesEqualTargetEntriesInTheWrittenDatabase()
    {
        // The assertion that matters, and the one no earlier test made: count what survives the
        // WRITE. Everything the suite asserted about entrapment proteins was true of the in-memory
        // list and still produced a database with 2.56 entrapment entries per target, because
        // ProteinDbWriter persists one entry per consensus and entrapment proteins are each their
        // own consensus.
        List<Protein> targets = ExpandedVariantProteins();
        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateEntrapment(
            targets, Tryptic, NothingForbidden);

        string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
            "entrapment_entry_parity.xml");
        ProteinDbWriter.WriteXmlDatabase(null, targets.Concat(entrapment).ToList(), path);

        string written = File.ReadAllText(path);
        int entries = System.Text.RegularExpressions.Regex.Matches(written, "<entry ").Count;
        int entrapmentEntries = System.Text.RegularExpressions.Regex
            .Matches(written, "<accession>Random_").Count;

        Assert.That(entries, Is.EqualTo(2), "one target entry and one entrapment entry");
        Assert.That(entrapmentEntries, Is.EqualTo(1));
        Assert.That(entries - entrapmentEntries, Is.EqualTo(entrapmentEntries),
            "the database must hold one entrapment entry per target entry");

        File.Delete(path);
    }

    [Test]
    public void EveryFoldGetsOneEntryPerTarget()
    {
        List<Protein> targets = ExpandedVariantProteins();

        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateEntrapment(
            targets, Tryptic, NothingForbidden, foldCount: 3);

        Assert.That(entrapment.Count, Is.EqualTo(3));
        Assert.That(entrapment.Select(p => p.Accession).ToList(),
            Is.EquivalentTo(new[] { "Random_P00001_f0", "Random_P00001_f1", "Random_P00001_f2" }));
    }

    [Test]
    public void AFastaStyleListWithNoVariantsIsUnchangedByTheEntryRule()
    {
        // The FASTA path was always correct, and must stay that way: no consensus sharing, so every
        // protein is its own entry and every one gets a partner.
        var targets = new List<Protein>
        {
            new Protein(Sequence, "P00001"),
            new Protein(Sequence, "P00002"),
        };

        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateEntrapment(
            targets, Tryptic, NothingForbidden);

        Assert.That(entrapment.Count, Is.EqualTo(2));
        Assert.That(entrapment.Select(p => p.Accession).ToList(),
            Is.EquivalentTo(new[] { "Random_P00001_f0", "Random_P00002_f0" }));
    }

    private const string ForeignSequence = "MQVLGYTTPDNRAWEDSFLGKQPTMLNDVAHER";

    [Test]
    public void GenerateEntrapmentCanReportEveryAssemblyAsItGoes()
    {
        // One loop rather than two. A caller building a QC report used to walk DatabaseEntries and
        // call Create itself, and those two loops had to agree about folds and about what an entry
        // is -- until one compared its partner count against the loaded protein count instead of the
        // entry count and failed silently for a day.
        var targets = new List<Protein>
        {
            new Protein(Sequence, "P00001"),
            new Protein("M" + new string('Q', 79), "Q156A1"),   // contributes no partner at all
        };

        var seen = new List<(string Accession, int Fold, bool Produced)>();
        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateEntrapment(
            targets, Tryptic, NothingForbidden,
            (entry, fold, assembly) => seen.Add((entry.Accession, fold, assembly.EntrapmentSequence.Length > 0)),
            foldCount: 2);

        Assert.That(seen.Select(x => x.Accession), Is.EquivalentTo(
            new[] { "P00001", "P00001", "Q156A1", "Q156A1" }), "every entry, every fold");
        Assert.That(seen.Where(x => x.Accession == "P00001").Select(x => x.Fold),
            Is.EquivalentTo(new[] { 0, 1 }));
        Assert.That(seen.Single(x => x.Accession == "Q156A1" && x.Fold == 0).Produced, Is.False,
            "an entry that contributes nothing is still reported -- that is what a ratio needs");
        Assert.That(entrapment, Has.Count.EqualTo(2), "and it is still absent from the output");
    }

    [Test]
    public void ANegativeSeedGivesTheSameDatabaseWhateverTheCulture()
    {
        // Interpolating the seed formats it through the current culture, and a negative one renders
        // its sign as U+002D under en-US but U+2212 MINUS SIGN under sv-SE, fi-FI and lt-LT --
        // different bytes into SHA-256, so the same request produced a different database on a
        // differently-configured machine. The reproducibility this construction exists for has to
        // survive a culture change, not only a process restart.
        CultureInfo original = CultureInfo.CurrentCulture;
        try
        {
            CultureInfo.CurrentCulture = new CultureInfo("en-US");
            string invariant = EntrapmentAssembler
                .Assemble(Sequence, Tryptic, NothingForbidden, seed: -7).EntrapmentSequence;

            foreach (string culture in new[] { "sv-SE", "fi-FI", "lt-LT" })
            {
                CultureInfo.CurrentCulture = new CultureInfo(culture);
                Assert.That(EntrapmentAssembler.Assemble(Sequence, Tryptic, NothingForbidden, seed: -7)
                    .EntrapmentSequence, Is.EqualTo(invariant), culture);
            }
        }
        finally
        {
            CultureInfo.CurrentCulture = original;
        }
    }

    [Test]
    public void AnAgentWhoseSitesCannotBeHeldIsRefused()
    {
        // trypsin|P's K[P]| pins nothing: CleavageSitePositions skips a position where a motif fits
        // but is prevented, so a rearrangement can move that P away and invent a cleavage site the
        // target lacked. Five shipped agents carry preventing motifs. Refusing is the honest answer
        // until pinning works across piece boundaries -- a database that digests differently from
        // its target produces peptides that fail to pair, and it does so silently.
        var preventing = new DigestionParams("trypsin|P", minPeptideLength: 7, maxMissedCleavages: 2);

        var ex = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentAssembler.Assemble(Sequence, preventing, NothingForbidden));
        Assert.That(ex.Message, Does.Contain("preventing-cleavage motif"));
    }

    [Test]
    public void AMultiResidueAgentIsRefused()
    {
        // StcE's TX|T spans three residues and is cut by the piece boundary, so one piece ends "TX"
        // and the next begins "T" -- neither fragment matches inside its own piece and nothing is
        // held at the seam.
        var multiResidue = new DigestionParams("StcE-trypsin", minPeptideLength: 7, maxMissedCleavages: 2);

        var ex = Assert.Throws<MzLibUtil.MzLibException>(() =>
            EntrapmentAssembler.Assemble(Sequence, multiResidue, NothingForbidden));
        Assert.That(ex.Message, Does.Contain("multi-residue motif"));
    }

    [Test]
    public void AProteinWhoseEveryPieceIsExcisedIsNotEmitted()
    {
        // Q156A1 is the real case: a methionine followed by 79 glutamines. Its single base piece has
        // exactly one arrangement once the termini are anchored, so it is excised and nothing
        // remains. Writing an empty-sequence protein makes a loader warn and discard it -- and its
        // decoy with it -- leaving a database whose entrapment count is quietly short of its target
        // count, which is exactly how the off-by-one at entry level was noticed.
        var homopolymer = new Protein("M" + new string('Q', 79), "Q156A1");
        var ordinary = new Protein(Sequence, "P00001");

        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateEntrapment(
            new[] { homopolymer, ordinary }, Tryptic, NothingForbidden);

        Assert.That(entrapment.Select(p => p.Accession), Is.EquivalentTo(new[] { "Random_P00001_f0" }),
            "the homopolymer contributes no entry at all, rather than an empty one");
        Assert.That(entrapment.All(p => p.BaseSequence.Length > 0), Is.True);
    }

    [Test]
    public void TheHomopolymerReallyHasNoArrangement()
    {
        // Guard on the premise of the test above: if anchoring ever stopped freezing this sequence,
        // that test would pass for the wrong reason.
        EntrapmentAssembly assembly = EntrapmentAssembler.Assemble("M" + new string('Q', 79),
            Tryptic, NothingForbidden);

        Assert.That(assembly.EntrapmentSequence, Is.Empty);
        Assert.That(assembly.ExcisedCount, Is.EqualTo(1));
        Assert.That(assembly.Pieces.Single().Failure,
            Is.EqualTo(EntrapmentFailure.NoPermutationExists));
    }

    [Test]
    public void AForeignAccessionRoundTripsAndNamesNoTarget()
    {
        string accession = EntrapmentAccession.FormatForeign("Q9XYZ1");

        Assert.That(accession, Is.EqualTo("Random_foreign_Q9XYZ1"));
        Assert.That(EntrapmentAccession.TryParseForeign(accession, out string foreign), Is.True);
        Assert.That(foreign, Is.EqualTo("Q9XYZ1"));

        // It is entrapment, so a consumer scanning for the identifier must still find it -- missing
        // one would count a foreign entry as a TARGET, which corrupts the estimate in the worst
        // direction.
        Assert.That(accession, Does.StartWith("Random_"));

        // But there is no target to name, so the target parser must refuse rather than guess.
        Assert.That(EntrapmentAccession.TryParse(accession, out _, out _), Is.False);
    }

    [Test]
    public void TargetParsingRefusesAForeignAccessionThatLooksLikeAFold()
    {
        // The edge case the marker exists for: a foreign protein whose own accession ends in
        // something shaped like a fold suffix. Refusing on the marker rather than the shape keeps
        // "foreign_ABC" from ever being reported as a target accession.
        string accession = EntrapmentAccession.FormatForeign("ABC_f1");

        Assert.That(accession, Is.EqualTo("Random_foreign_ABC_f1"));
        Assert.That(EntrapmentAccession.TryParse(accession, out string target, out _), Is.False,
            "a foreign entry has no target, whatever its accession happens to end with");
        Assert.That(target, Is.Empty);
        Assert.That(EntrapmentAccession.TryParseForeign(accession, out string foreign), Is.True);
        Assert.That(foreign, Is.EqualTo("ABC_f1"));
    }

    [Test]
    public void AForeignEntryIsRelabelledNotRearranged()
    {
        // A foreign protein already is a sequence the sample cannot contain, so permuting it would
        // buy nothing and lose the annotations that describe it.
        var foreign = new Protein(ForeignSequence, "Q9XYZ1",
            sequenceVariations: new List<SequenceVariation>
            {
                new SequenceVariation(5, 5, "G", "A", "irrelevant to entrapment"),
            });

        Protein entry = EntrapmentProteinGenerator.CreateForeign(foreign);

        Assert.That(entry.BaseSequence, Is.EqualTo(ForeignSequence), "the sequence is not touched");
        Assert.That(entry.IsEntrapment, Is.True);
        Assert.That(entry.Accession, Is.EqualTo("Random_foreign_Q9XYZ1"));
        Assert.That(entry.SequenceVariations, Is.Empty,
            "applying variants would expand one entry into several, and the arm is one entry per entry");
    }

    [Test]
    public void ForeignPeptidesSharedWithTheTargetAreReported()
    {
        // Homology is the hazard: a conserved protein shares peptides with its ortholog, and a
        // shared peptide is a REAL target peptide sitting in the entrapment database. Nothing can be
        // permuted away here, so it is counted.
        IDigestionParams digestion = Tryptic;
        var target = new Protein(Sequence, "P00001");

        var targetPeptides = target.Digest(digestion, new List<Modification>(), new List<Modification>())
            .Select(pep => pep.BaseSequence).ToHashSet();
        string conserved = targetPeptides.First(pep => pep.Length >= 9);

        // one foreign protein that shares a peptide, and one that does not
        var foreign = new List<Protein>
        {
            new Protein("MQVLGYTTPDNR" + conserved + "AWEDSFLGK", "Q9SHARED"),
            new Protein(ForeignSequence, "Q9CLEAN"),
        };

        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateForeignEntrapment(
            foreign, digestion, targetPeptides, out var shared);

        Assert.That(entrapment.Select(e => e.Accession),
            Is.EquivalentTo(new[] { "Random_foreign_Q9SHARED", "Random_foreign_Q9CLEAN" }));
        Assert.That(shared.ContainsKey("Q9SHARED"), Is.True,
            "fixture must actually share a peptide, or it proves nothing");
        Assert.That(shared["Q9SHARED"], Does.Contain(conserved));
        Assert.That(shared.ContainsKey("Q9CLEAN"), Is.False);
    }

    [Test]
    public void DisjointForeignProteomeSharesNothing()
    {
        IDigestionParams digestion = Tryptic;
        var target = new Protein(Sequence, "P00001");
        var targetPeptides = target.Digest(digestion, new List<Modification>(), new List<Modification>())
            .Select(pep => pep.BaseSequence).ToHashSet();

        List<Protein> entrapment = EntrapmentProteinGenerator.GenerateForeignEntrapment(
            new[] { new Protein(ForeignSequence, "Q9CLEAN") }, digestion, targetPeptides, out var shared);

        Assert.That(entrapment, Has.Count.EqualTo(1));
        Assert.That(shared, Is.Empty);
    }

    private static string Sorted(string s) => string.Concat(s.OrderBy(c => c));

    // Two cysteines, well apart, so a disulfide bond has somewhere to be mapped to.
    private const string ProteoformSequence = "MADCQVLGYTTPDNRAWEDSFLCGKQPTMLNDVAHERGGLYT";

    [Test]
    public void AProteoformPartnerIsTheWholeProteinRearranged()
    {
        var target = new Protein(ProteoformSequence, "P00001");

        Protein partner = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden,
            out EntrapmentAssembly assembly);

        Assert.That(partner, Is.Not.Null);
        Assert.That(partner.BaseSequence, Is.Not.EqualTo(ProteoformSequence));
        Assert.That(Sorted(partner.BaseSequence), Is.EqualTo(Sorted(ProteoformSequence)),
            "isomeric with its target, as everywhere else");
        Assert.That(partner.BaseSequence, Has.Length.EqualTo(ProteoformSequence.Length),
            "nothing is excised in this mode, which is what makes annotations carryable");
        Assert.That(partner.IsEntrapment, Is.True);
        Assert.That(assembly.ExcisedCount, Is.Zero);

        // The protein's own termini stay put, so terminal modifications survive.
        Assert.That(partner.BaseSequence[0], Is.EqualTo(ProteoformSequence[0]));
        Assert.That(partner.BaseSequence[^1], Is.EqualTo(ProteoformSequence[^1]));
    }

    [Test]
    public void ProteoformModeIgnoresCleavageSites()
    {
        // The difference from the bottom-up path: top-down does not digest, so K and R are ordinary
        // residues and are free to move. If they were still pinned this would be the same thing
        // twice under two names.
        var target = new Protein(ProteoformSequence, "P00001");

        Protein partner = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden, out _);

        var targetKr = Enumerable.Range(0, ProteoformSequence.Length)
            .Where(i => ProteoformSequence[i] is 'K' or 'R').ToList();
        var partnerKr = Enumerable.Range(0, partner.BaseSequence.Length)
            .Where(i => partner.BaseSequence[i] is 'K' or 'R').ToList();

        Assert.That(partnerKr, Is.Not.EqualTo(targetKr),
            "cleavage residues are not held in proteoform mode");
    }

    [Test]
    public void ADisulfideBondStillJoinsTheSameTwoCysteines()
    {
        // The reason annotations are mapped rather than copied. Copying the coordinates would assert
        // a bond between whatever residues happen to land there, which is simply false.
        int first = ProteoformSequence.IndexOf('C') + 1;
        int second = ProteoformSequence.LastIndexOf('C') + 1;
        Assert.That(first, Is.LessThan(second), "fixture needs two cysteines");

        var target = new Protein(ProteoformSequence, "P00001",
            disulfideBonds: new List<DisulfideBond> { new(first, second, "test bond") });

        Protein partner = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden, out _);

        Assert.That(partner.DisulfideBonds, Has.Count.EqualTo(1));
        DisulfideBond bond = partner.DisulfideBonds.Single();
        Assert.That(partner.BaseSequence[bond.OneBasedBeginPosition - 1], Is.EqualTo('C'));
        Assert.That(partner.BaseSequence[bond.OneBasedEndPosition - 1], Is.EqualTo('C'));
        Assert.That(bond.OneBasedBeginPosition, Is.LessThan(bond.OneBasedEndPosition),
            "a rearrangement can put the later residue first, so the order must be restored");
    }

    [Test]
    public void TruncationProductsKeepTheirSpans()
    {
        // A span has no contiguous image under a rearrangement, and the span is what defines a
        // proteoform: the entrapment proteoform truncated at 25 is the counterpart of the target
        // truncated at 25. Its residues differ, which is the entire point of an entrapment sequence.
        var target = new Protein(ProteoformSequence, "P00001",
            proteolysisProducts: new List<TruncationProduct> { new(2, 25, "chain") });

        Protein partner = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden, out _);

        Assert.That(partner.TruncationProducts, Has.Count.EqualTo(1));
        TruncationProduct product = partner.TruncationProducts.Single();
        Assert.That(product.OneBasedBeginPosition, Is.EqualTo(2));
        Assert.That(product.OneBasedEndPosition, Is.EqualTo(25));
        Assert.That(product.OneBasedEndPosition, Is.LessThanOrEqualTo(partner.BaseSequence.Length),
            "the coordinate must stay in range -- out-of-range features are what made an earlier "
            + "generation of these databases unloadable");
    }

    [Test]
    public void ProteoformModifsAreCarriedNotLost()
    {
        var mods = new Dictionary<int, List<Modification>>
        {
            // Phospho targets T, so it must sit on one -- a modification whose motif does not
            // fit its position is invalid and mzLib drops it, which is correct and was what
            // this fixture originally proved rather than what it meant to test.
            { 10, new List<Modification> { Phospho() } },
        };
        var target = new Protein(ProteoformSequence, "P00001", oneBasedModifications: mods);

        Protein partner = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden, out _);

        int carried = partner.OneBasedPossibleLocalizedModifications.Sum(kv => kv.Value.Count);
        Assert.That(carried, Is.EqualTo(1), "nothing is excised, so no modification has anywhere to fall");
    }

    [Test]
    public void AProteoformWithNoRearrangementIsDroppedWhole()
    {
        // Half a proteoform is not a proteoform. Where the bottom-up path excises one piece and
        // keeps the rest, this returns nothing at all.
        var target = new Protein("AAAAAAAA", "P00001");

        Protein partner = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden,
            out EntrapmentAssembly assembly);

        Assert.That(partner, Is.Null);
        Assert.That(assembly.ExcisedCount, Is.EqualTo(1));
    }

    [Test]
    public void ProteoformModeIsDeterministic()
    {
        var target = new Protein(ProteoformSequence, "P00001");

        Protein first = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden, out _);
        Protein second = EntrapmentProteinGenerator.CreateProteoform(target, NothingForbidden, out _);

        Assert.That(second.BaseSequence, Is.EqualTo(first.BaseSequence));
    }
}
