using System;
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
        foreach (TruncationProduct tp in entrapment.TruncationProducts)
        {
            Assert.That(tp.OneBasedEndPosition, Is.LessThanOrEqualTo(entrapment.BaseSequence.Length));
        }
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
}
