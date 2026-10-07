using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using Chemistry;
using Omics;
using Omics.SequenceConversion;
using Readers.ProForma;
using System;
using System.Linq;

namespace Test.Omics.SequenceConversion;
[TestFixture]
public class PeptideAndProteinConversion
{
    #region PeptideWithSetModifications Conversion Tests

    [Test]
    public static void TestConvertModificationsOnPeptideWithSetModifications()
    {
        // Create a protein with MetaMorpheus-style modifications
        ModificationMotif.TryGetMotif("S", out var motifS);
        ModificationMotif.TryGetMotif("K", out var motifK);

        var phosphoMM = new Modification(
            _originalId: "Phosphorylation",
            _modificationType: "Common Biological",
            _target: motifS,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

        var acetylMM = new Modification(
            _originalId: "Acetylation",
            _modificationType: "Common Biological",
            _target: motifK,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"));

        var protein = new Protein("PEPTKSDE", "TestProtein");

        var modsOneIsNterm = new Dictionary<int, Modification>
           {
               { 1, acetylMM }, // N-terminal acetylation
               { 4, acetylMM }, // Acetyl on K at position 4 (P-E-P-T-K)
               { 6, phosphoMM }  // Phospho on S at position 6 (P-E-P-T-K-S)
           };

        var digestionParams = new DigestionParams(protease: "trypsin");
        var peptide = new PeptideWithSetModifications(
            protein,
            digestionParams,
            oneBasedStartResidueInProtein: 1,
            oneBasedEndResidueInProtein: 8,
            cleavageSpecificity: CleavageSpecificity.Full,
            peptideDescription: "Test",
            missedCleavages: 0,
            allModsOneIsNterminus: modsOneIsNterm,
            numFixedMods: 0);

        // Verify original mods are MetaMorpheus style
        Assert.That(peptide.AllModsOneIsNterminus[1].ModificationType, Is.EqualTo("Common Biological"));
        Assert.That(peptide.AllModsOneIsNterminus[4].ModificationType, Is.EqualTo("Common Biological"));
        Assert.That(peptide.AllModsOneIsNterminus[6].ModificationType, Is.EqualTo("Common Biological"));

        // Convert to Unimod convention
        peptide.ConvertModifications(UnimodModificationLookup.Instance);

        // Verify conversions
        Assert.That(peptide.AllModsOneIsNterminus[1].ModificationType, Is.EqualTo("Unimod"));
        Assert.That(peptide.AllModsOneIsNterminus[4].ModificationType, Is.EqualTo("Unimod"));
        Assert.That(peptide.AllModsOneIsNterminus[6].ModificationType, Is.EqualTo("Unimod"));

        // Verify chemical formulas are preserved
        Assert.That(peptide.AllModsOneIsNterminus[1].ChemicalFormula, Is.Not.Null);
        Assert.That(peptide.AllModsOneIsNterminus[4].ChemicalFormula, Is.Not.Null);
        Assert.That(peptide.AllModsOneIsNterminus[6].ChemicalFormula, Is.Not.Null);
    }

    [Test]
    public static void TestConvertModificationsOnPeptidePreservesMotifs()
    {
        // Test that conversion preserves amino acid targets
        ModificationMotif.TryGetMotif("M", out var motifM);

        var oxidationMM = new Modification(
            _originalId: "Oxidation",
            _modificationType: "Common Variable",
            _target: motifM,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("O1"));

        var protein = new Protein("PEPTMIDE", "TestProtein");

        var modsOneIsNterm = new Dictionary<int, Modification>
           {
               { 7, oxidationMM } // Oxidation on M
           };

        var digestionParams = new DigestionParams(protease: "trypsin");
        var peptide = new PeptideWithSetModifications(
            protein,
            digestionParams,
            oneBasedStartResidueInProtein: 1,
            oneBasedEndResidueInProtein: 8,
            cleavageSpecificity: CleavageSpecificity.Full,
            peptideDescription: "Test",
            missedCleavages: 0,
            allModsOneIsNterminus: modsOneIsNterm,
            numFixedMods: 0);

        // Get original target
        var originalTarget = peptide.AllModsOneIsNterminus[7].Target.ToString();

        // Convert to Unimod
        peptide.ConvertModifications(UnimodModificationLookup.Instance);

        // Verify target is preserved
        Assert.That(peptide.AllModsOneIsNterminus[7].Target.ToString(), Is.EqualTo(originalTarget));
        Assert.That(peptide.AllModsOneIsNterminus[7].Target.ToString(), Does.Contain("M"));
    }

    [Test]
    public static void TestConvertModificationsOnPeptideWithCTerminalMod()
    {
        // Test conversion with C-terminal modification
        ModificationMotif.TryGetMotif("X", out var motifX);

        var amidationMM = new Modification(
            _originalId: "Amidation",
            _modificationType: "Common Biological",
            _target: motifX,
            _locationRestriction: "Peptide C-terminal.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H1N1"));

        var protein = new Protein("PEPTIDE", "TestProtein");

        var modsOneIsNterm = new Dictionary<int, Modification>
           {
               { 9, amidationMM } // C-terminal mod (length 7 + 2)
           };

        var digestionParams = new DigestionParams(protease: "trypsin");
        var peptide = new PeptideWithSetModifications(
            protein,
            digestionParams,
            oneBasedStartResidueInProtein: 1,
            oneBasedEndResidueInProtein: 7,
            cleavageSpecificity: CleavageSpecificity.Full,
            peptideDescription: "Test",
            missedCleavages: 0,
            allModsOneIsNterminus: modsOneIsNterm,
            numFixedMods: 0);

        // Convert to Unimod
        peptide.ConvertModifications(UnimodModificationLookup.Instance);

        // Verify C-terminal mod was converted
        Assert.That(peptide.AllModsOneIsNterminus.ContainsKey(9), Is.True);
        Assert.That(peptide.AllModsOneIsNterminus[9].ModificationType, Is.EqualTo("Unimod"));
    }

    [Test]
    public static void TestConvertModificationsOnEmptyPeptide()
    {
        // Test that conversion works with no modifications
        var protein = new Protein("PEPTIDE", "TestProtein");

        var modsOneIsNterm = new Dictionary<int, Modification>();

        var digestionParams = new DigestionParams(protease: "trypsin");
        var peptide = new PeptideWithSetModifications(
            protein,
            digestionParams,
            oneBasedStartResidueInProtein: 1,
            oneBasedEndResidueInProtein: 7,
            cleavageSpecificity: CleavageSpecificity.Full,
            peptideDescription: "Test",
            missedCleavages: 0,
            allModsOneIsNterminus: modsOneIsNterm,
            numFixedMods: 0);

        // Should not throw
        Assert.DoesNotThrow(() => peptide.ConvertModifications(UnimodModificationLookup.Instance));

        // Should still have no mods
        Assert.That(peptide.AllModsOneIsNterminus.Count, Is.EqualTo(0));
    }

    [Test]
    public static void TestPeptideConversionRoundTrip()
    {
        // Test converting from MetaMorpheus to Unimod and back preserves chemistry
        ModificationMotif.TryGetMotif("S", out var motifS);

        var phosphoMM = new Modification(
            _originalId: "Phosphorylation",
            _modificationType: "Common Biological",
            _target: motifS,
            _locationRestriction: "Anywhere.",
            _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));

        var protein = new Protein("PEPTSIDE", "TestProtein");

        var modsOneIsNterm = new Dictionary<int, Modification>
           {
               { 7, phosphoMM }
           };

        var digestionParams = new DigestionParams(protease: "trypsin");
        var peptide = new PeptideWithSetModifications(
            protein,
            digestionParams,
            oneBasedStartResidueInProtein: 1,
            oneBasedEndResidueInProtein: 8,
            cleavageSpecificity: CleavageSpecificity.Full,
            peptideDescription: "Test",
            missedCleavages: 0,
            allModsOneIsNterminus: modsOneIsNterm,
            numFixedMods: 0);

        var originalFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;
        var originalTarget = peptide.AllModsOneIsNterminus[7].Target.ToString();

        // Convert to Unimod
        peptide.ConvertModifications(UnimodModificationLookup.Instance);

        var UnimodFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;

        // Convert back to MetaMorpheus
        peptide.ConvertModifications(MzLibModificationLookup.Instance);

        var finalFormula = peptide.AllModsOneIsNterminus[7].ChemicalFormula;
        var finalTarget = peptide.AllModsOneIsNterminus[7].Target.ToString();

        // Verify chemistry is preserved
        Assert.That(originalFormula.Equals(UnimodFormula), Is.True);
        Assert.That(originalFormula.Equals(finalFormula), Is.True);
        Assert.That(originalTarget, Is.EqualTo(finalTarget));
    }

    #endregion

    #region Full sequences written for digested peptides

    /// <summary>
    /// Peptides digested from proteins that carry each protein modification mzLib's catalogs define (MetaMorpheus,
    /// TMT, UniProt, UNIMOD) where its motif fits: as a variable modification, and as a localized one on a target
    /// and on a decoy protein (decoys put localized modifications on residues whatever their location restriction).
    /// </summary>
    private static readonly Lazy<List<PeptideWithSetModifications>> DigestedCatalogPeptides = new(() =>
    {
        var peptides = new List<PeptideWithSetModifications>();
        var digestionParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1,
            maxModsForPeptides: 2, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
        var catalogs = Mods.MetaMorpheusProteinModifications.Concat(Mods.IsobaricLabelModifications)
            .Concat(Mods.UniprotModifications).Concat(Mods.UnimodModifications);
        foreach (var mod in catalogs.Where(m => m.Target != null && m.MonoisotopicMass.HasValue && m.ValidModification))
        {
            var residues = new string(mod.Target.ToString().Select(c => c is 'X' or 'x' ? 'A' : char.ToUpperInvariant(c)).ToArray());
            var sequence = residues + "GGGGK" + residues + "GGGGK" + "GGG" + residues;
            Dictionary<int, List<Modification>> Localized() => new[] { 1, residues.Length + 6, sequence.Length }
                .ToDictionary(position => position, _ => new List<Modification> { mod });
            var proteins = new[]
            {
                (Protein: new Protein(sequence, "P"), Variable: new List<Modification> { mod }),
                (Protein: new Protein(sequence, "P", oneBasedModifications: Localized()), Variable: new List<Modification>()),
                (Protein: new Protein(sequence, "DECOY_P", isDecoy: true, oneBasedModifications: Localized()), Variable: new List<Modification>())
            };
            foreach (var (protein, variable) in proteins)
            {
                peptides.AddRange(protein.Digest(digestionParams, new List<Modification>(), variable)
                    .Cast<PeptideWithSetModifications>()
                    .Where(p => p.AllModsOneIsNterminus.Count > 0));
            }
        }
        return peptides;
    });

    private static int KeyOf(CanonicalModification mod, int length) => mod.PositionType switch
    {
        ModificationPositionType.NTerminus => 1,
        ModificationPositionType.CTerminus => length + 2,
        _ => mod.ResidueIndex!.Value + 2
    };

    [Test]
    public static void MzLibParser_DigestedCatalogModification_CarriesItsCatalogEntryAndUnimodId()
    {
        int identical = 0, sameNamedEntry = 0, withoutMismatchedId = 0;
        foreach (var peptide in DigestedCatalogPeptides.Value)
        {
            var parsed = MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value;
            Assert.That(parsed.Modifications.Length, Is.EqualTo(peptide.AllModsOneIsNterminus.Count), peptide.FullSequence);
            foreach (var mod in parsed.Modifications)
            {
                var original = peptide.AllModsOneIsNterminus[KeyOf(mod, parsed.BaseSequence.Length)];
                var attached = mod.MzLibModification;
                Assert.That(attached, Is.Not.Null, peptide.FullSequence);
                var cited = CanonicalModification.GetUnimodId(original);
                if (cited.HasValue && mod.UnimodId == null)
                {
                    // Left out only when the cited UNIMOD record is another chemical (DiLeu-12plex citing 1327, UniProt
                    // N,N-dimethylproline citing 529, ...).
                    var recordMass = Mods.UnimodModifications.First(m => m.ModificationType == "Unimod" && CanonicalModification.GetUnimodId(m) == cited).MonoisotopicMass!.Value;
                    Assert.That(Math.Abs(recordMass - original.MonoisotopicMass!.Value), Is.GreaterThan(1), peptide.FullSequence);
                    withoutMismatchedId++;
                }
                else
                    Assert.That(mod.UnimodId, Is.EqualTo(cited), peptide.FullSequence);
                if (ReferenceEquals(attached, original))
                {
                    identical++;
                    continue;
                }

                // The name can't tell apart entries that share it and differ only in location restriction
                // (UNIMOD's N-terminal and peptide N-terminal Acetyl on X); the entry still fits where it was read.
                sameNamedEntry++;
                Assert.That(attached!.ModificationType, Is.EqualTo(original.ModificationType), peptide.FullSequence);
                Assert.That(attached.IdWithMotif, Is.EqualTo(original.IdWithMotif), peptide.FullSequence);
                Assert.That(attached.MonoisotopicMass, Is.EqualTo(original.MonoisotopicMass), peptide.FullSequence);
                // On a residue, digestion leaves terminal modifications of decoys (and protease products) where they
                // were: an N-terminal one on the first residue, a C-terminal one on the last, keeps that class unless the
                // name is also an entry allowed on any residue, which is what a residue is then read as.
                var unrestricted = mod.PositionType == ModificationPositionType.Residue
                                   && !ModificationLocalization.IsNTerminal(attached) && !ModificationLocalization.IsCTerminal(attached);
                var firstResidueNTerminal = mod.PositionType == ModificationPositionType.Residue && mod.ResidueIndex == 0
                                            && ModificationLocalization.IsNTerminal(original);
                var lastResidueCTerminal = mod.PositionType == ModificationPositionType.Residue && mod.ResidueIndex == parsed.BaseSequence.Length - 1
                                           && ModificationLocalization.IsCTerminal(original);
                if (mod.PositionType == ModificationPositionType.NTerminus || firstResidueNTerminal)
                    Assert.That(ModificationLocalization.IsNTerminal(attached) || unrestricted, peptide.FullSequence);
                if (mod.PositionType == ModificationPositionType.CTerminus || lastResidueCTerminal)
                    Assert.That(ModificationLocalization.IsCTerminal(attached) || unrestricted, peptide.FullSequence);
            }
        }
        Assert.That(identical, Is.GreaterThan(0));
        Assert.That(sameNamedEntry, Is.LessThan(identical));
        Assert.That(withoutMismatchedId, Is.GreaterThan(0));
    }

    [TestCase("Less Common", "Methylation on X", "N-terminal.")]
    [TestCase("Less Common", "Methylation on X", "C-terminal.")]
    [TestCase("Unimod", "Methyl on X", "N-terminal.")]
    [TestCase("Unimod", "Methyl on X", "Peptide C-terminal.")]
    [TestCase("Unimod", "Ethyl on X", "N-terminal.")]
    [TestCase("Unimod", "Ethyl on X", "Peptide C-terminal.")]
    [TestCase("Unimod", "Propyl on X", "Peptide N-terminal.")]
    [TestCase("Unimod", "Propyl on X", "C-terminal.")]
    public static void MzLibParserAndBuilder_NameSharedByTerminalEntries_KeepTheEntryForItsTerminus(string type, string id, string restriction)
    {
        // Mods.txt defines "Less Common:Methylation on X" twice and UNIMOD has several "Methyl on X"; reading the
        // full sequence through a dictionary keyed by IdWithMotif gives every one of them the same entry.
        var mod = Mods.AllProteinModsList.Single(m => m.ModificationType == type && m.IdWithMotif == id && m.LocationRestriction == restriction);
        var digestionParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1);
        var peptide = new Protein("AGGGGKAGGGGKGGGA", "P")
            .Digest(digestionParams, new List<Modification>(), new List<Modification> { mod })
            .Cast<PeptideWithSetModifications>()
            .First(p => p.AllModsOneIsNterminus.Count == 1);
        var (key, original) = peptide.AllModsOneIsNterminus.Single();

        var parsed = MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value;
        var attached = parsed.Modifications.Single().MzLibModification!;
        var built = PeptideWithSetModifications.FromCanonicalSequence(parsed);

        Assert.That(ModificationLocalization.IsNTerminal(attached), Is.EqualTo(ModificationLocalization.IsNTerminal(original)));
        Assert.That(ModificationLocalization.IsCTerminal(attached), Is.EqualTo(ModificationLocalization.IsCTerminal(original)));
        Assert.That(built.FullSequence, Is.EqualTo(peptide.FullSequence));
        Assert.That(built.AllModsOneIsNterminus.Keys, Is.EqualTo(new[] { key }));
        Assert.That(built.AllModsOneIsNterminus[key], Is.SameAs(attached));
    }

    [Test]
    public static void FromCanonicalSequence_DigestedPeptides_MatchTheOriginalPeptide()
    {
        // Targets and decoys alike, including the C-terminal modifications decoys and protease products keep on the
        // last residue. The objects differ only where entries share a name (see the parser test).
        int digested = 0;
        foreach (var peptide in DigestedCatalogPeptides.Value)
        {
            var parsed = MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value;
            var built = PeptideWithSetModifications.FromCanonicalSequence(parsed);

            Assert.That(built.BaseSequence, Is.EqualTo(peptide.BaseSequence));
            Assert.That(built.FullSequence, Is.EqualTo(peptide.FullSequence));
            Assert.That(built.MonoisotopicMass, Is.EqualTo(peptide.MonoisotopicMass).Within(1e-6), peptide.FullSequence);
            Assert.That(built.AllModsOneIsNterminus.Keys.Order(), Is.EqualTo(peptide.AllModsOneIsNterminus.Keys.Order()), peptide.FullSequence);
            foreach (var (key, mod) in built.AllModsOneIsNterminus)
                Assert.That(mod, Is.SameAs(parsed.Modifications.Single(m => KeyOf(m, parsed.BaseSequence.Length) == key).MzLibModification));
            digested++;
        }
        Assert.That(digested, Is.GreaterThan(0));
    }

    [Test]
    public static void FromCanonicalSequence_DecoyWithCTerminalModificationsOnTheLastResidueAndTheCTerminus_KeepsBoth()
    {
        // A decoy keeps its localized C-terminal amide on the last residue while a variable C-terminal modification
        // takes the C-terminus; both stay where they were written.
        var amide = Mods.UniprotModifications.Single(m => m.IdWithMotif == "Arginine amide on R");
        var amidation = Mods.MetaMorpheusProteinModifications.Single(m => m.IdWithMotif == "Amidation on X");
        var localized = new Dictionary<int, List<Modification>> { [12] = new() { amide } };
        var peptide = new Protein("APEPTIDEKAAR", "DECOY_P", isDecoy: true, oneBasedModifications: localized)
            .Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 2),
                new List<Modification>(), new List<Modification> { amidation })
            .Cast<PeptideWithSetModifications>()
            .Single(p => p.AllModsOneIsNterminus.Count == 2);
        Assert.That(peptide.FullSequence, Is.EqualTo("AAR[UniProt:Arginine amide on R]-[Less Common:Amidation on X]"));

        var built = PeptideWithSetModifications.FromCanonicalSequence(MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value);

        Assert.That(built.FullSequence, Is.EqualTo(peptide.FullSequence));
        Assert.That(built.AllModsOneIsNterminus[4], Is.SameAs(amide));
        Assert.That(built.AllModsOneIsNterminus[5], Is.SameAs(amidation));
    }

    [Test]
    public static void FromCanonicalSequence_ProFormaWrittenPeptide_ResolvesThroughTheFallbackLookup()
    {
        // The ProForma parser carries UNIMOD ids, not modification objects, so every modification goes to a lookup.
        var oxidation = Mods.MetaMorpheusProteinModifications.Single(m => m.ModificationType == "Common Variable" && m.IdWithMotif == "Oxidation on M");
        var phospho = Mods.MetaMorpheusProteinModifications.Single(m => m.ModificationType == "Common Biological" && m.IdWithMotif == "Phosphorylation on S");
        var peptide = new Protein("PEPMSIDEK", "P")
            .Digest(new DigestionParams(protease: "trypsin", maxModsForPeptides: 2), new List<Modification>(), new List<Modification> { oxidation, phospho })
            .Cast<PeptideWithSetModifications>()
            .Single(p => p.AllModsOneIsNterminus.Count == 2);

        var proForma = peptide.ToProFormaString();
        var parsed = ProFormaSequenceParser.Instance.Parse(proForma)!.Value;
        Assert.That(parsed.Modifications.All(m => m.MzLibModification == null), proForma);

        var byDefault = PeptideWithSetModifications.FromCanonicalSequence(parsed);
        var byUnimod = PeptideWithSetModifications.FromCanonicalSequence(parsed, UnimodModificationLookup.Instance);

        foreach (var built in new[] { byDefault, byUnimod })
        {
            Assert.That(built.AllModsOneIsNterminus.Keys.Order(), Is.EqualTo(peptide.AllModsOneIsNterminus.Keys.Order()));
            Assert.That(built.MonoisotopicMass, Is.EqualTo(peptide.MonoisotopicMass).Within(1e-6));
        }
        Assert.That(byUnimod.AllModsOneIsNterminus.Values.Select(m => m.ModificationType), Is.All.EqualTo("Unimod"));
    }

    [Test]
    public static void ProFormaParser_DigestedCatalogModificationWrittenByAccession_CarriesItsCatalogEntry()
    {
        // Modifications without a usable UNIMOD id are written by PSI-MOD accession, from the peptide and from its
        // full sequence alike.
        int identical = 0, sameAccessionEntry = 0;
        foreach (var peptide in DigestedCatalogPeptides.Value)
        {
            var fromFullSequence = MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value;
            var written = new[]
            {
                (ProForma: peptide.ToProFormaString(), Mods: (IDictionary<int, Modification>)peptide.AllModsOneIsNterminus),
                (ProForma: ProFormaSequenceSerializer.Instance.Serialize(fromFullSequence),
                    Mods: fromFullSequence.Modifications.ToDictionary(m => KeyOf(m, fromFullSequence.BaseSequence.Length), m => m.MzLibModification!))
            };
            foreach (var (proForma, mods) in written.Where(w => w.ProForma!.Contains("[MOD:") || w.ProForma.Contains("[RESID:")))
            {
                var parsed = ProFormaSequenceParser.Instance.Parse(proForma!)!.Value;
                Assert.That(parsed.Modifications.Length, Is.EqualTo(mods.Count), proForma);
                foreach (var mod in parsed.Modifications)
                {
                    var original = mods[KeyOf(mod, parsed.BaseSequence.Length)];
                    var attached = mod.MzLibModification;
                    Assert.That(attached, Is.Not.Null, proForma);
                    Assert.That(mod.UnimodId, Is.Null, proForma);
                    if (ReferenceEquals(attached, original))
                    {
                        identical++;
                        continue;
                    }

                    // The accession can't tell apart entries that share it on one residue (UniProt lists MOD:00165 for
                    // N-linked (Hex) and N-linked (Man) tryptophan).
                    sameAccessionEntry++;
                    Assert.That(ProFormaConverter.BuildDescriptor(attached!).Value, Is.EqualTo(mod.OriginalRepresentation), proForma);
                    Assert.That(attached.Target.ToString(), Is.EqualTo(original.Target.ToString()), proForma);
                    Assert.That(attached.LocationRestriction, Is.EqualTo(original.LocationRestriction), proForma);
                    Assert.That(attached.MonoisotopicMass, Is.EqualTo(original.MonoisotopicMass), proForma);
                }
            }
        }
        Assert.That(identical, Is.GreaterThan(0));
        Assert.That(sameAccessionEntry, Is.LessThan(identical));
    }

    [Test]
    public static void FromCanonicalSequence_ModificationInNoCatalog_ThrowsNamingIt()
    {
        // A custom modification, read by MetaMorpheus from a user's file, is in none of mzLib's modification catalogs.
        ModificationMotif.TryGetMotif("K", out var motifK);
        var custom = new Modification(_originalId: "Nameless", _modificationType: "Custom", _target: motifK,
            _locationRestriction: "Anywhere.", _monoisotopicMass: 100.0);
        var peptide = new Protein("PEPKR", "P")
            .Digest(new DigestionParams(protease: "top-down", minPeptideLength: 1), new List<Modification>(), new List<Modification> { custom })
            .Cast<PeptideWithSetModifications>()
            .First(p => p.AllModsOneIsNterminus.Count > 0);
        var parsed = MzLibSequenceParser.Instance.Parse(peptide.FullSequence)!.Value;

        Assert.That(parsed.Modifications.Single().MzLibModification, Is.Null);
        Assert.That(() => PeptideWithSetModifications.FromCanonicalSequence(parsed),
            Throws.TypeOf<SequenceConversionException>().With.Message.Contains("Custom:Nameless on K"));
    }

    [Test]
    public static void GlobalModificationLookupProteinOnly_LeavesOutRnaModifications()
    {
        var rnaOnly = Mods.MetaMorpheusRnaModifications.First(m => Mods.AllProteinModsList.All(p => p.IdWithMotif != m.IdWithMotif));
        var rnaName = $"{rnaOnly.ModificationType}:{rnaOnly.IdWithMotif}";
        var rnaMod = CanonicalModification.AtResidue(0, rnaOnly.Target.ToString()[0], rnaName, mzLibId: rnaName);
        var acetyllysine = Mods.UniprotModifications.Single(m => m.IdWithMotif == "N6-acetyllysine on K");
        var proteinMod = CanonicalModification.AtResidue(0, 'K', "UniProt:N6-acetyllysine on K", mzLibId: "UniProt:N6-acetyllysine on K");

        Assert.That(GlobalModificationLookup.Instance.TryResolve(rnaMod)?.MzLibModification, Is.SameAs(rnaOnly));
        Assert.That(GlobalModificationLookup.ProteinOnly.TryResolve(rnaMod), Is.Null);
        Assert.That(GlobalModificationLookup.ProteinOnly.TryResolve(proteinMod)?.MzLibModification, Is.SameAs(acetyllysine));
    }

    #endregion
}
