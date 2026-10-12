using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Omics.SequenceConversion;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.Omics.SequenceConversion;

/// <summary>
/// Terminal modifications identified by UNIMOD accession, resolved against the real Unimod catalog.
/// This is what a ProForma "[UNIMOD:1]-PEPTIDE" or an mzIdentML Modification at location 0 hands the lookup:
/// the accession, the position, usually the residue it sits on, and often the mass delta.
///
/// Each case was first run against the unfixed lookup, which returned a modification whose Unimod id
/// merely contained the requested one, or a side-chain modification in place of the terminal one:
/// N-terminal Acetyl on M became "Met->Hse on M", on A "Fluoro on A", with no residue "Acetyl on T";
/// C-terminal Amidated on K became "Phospho on K".
///
/// A new lookup is built per case because the lookup caches per instance, null results included.
/// </summary>
[TestFixture]
public class UnimodModificationLookupTests
{
    [TestCase(1, null, null, "Acetyl on X")]
    [TestCase(1, null, 42.010565, "Acetyl on X")]
    [TestCase(1, 'M', null, "Acetyl on X")]
    [TestCase(1, 'M', 42.010565, "Acetyl on X")]
    [TestCase(1, 'A', 42.010565, "Acetyl on X")]
    [TestCase(1, 'S', null, "Acetyl on X")]
    [TestCase(1, 'K', 42.010565, "Acetyl on X")]
    [TestCase(4, 'C', 57.021464, "Carbamidomethyl on X")]
    [TestCase(5, 'K', null, "Carbamyl on X")]
    [TestCase(7, 'F', null, "Deamidated on F")]
    [TestCase(27, 'E', -18.010565, "Glu->pyro-Glu on E")]
    [TestCase(28, 'Q', null, "Gln->pyro-Glu on Q")]
    [TestCase(34, 'A', null, "Methyl on X")]
    [TestCase(36, 'P', null, "Dimethyl on P")]
    [TestCase(36, 'A', null, "Dimethyl on X")]
    [TestCase(122, 'G', null, "Formyl on X")]
    // Unimod 214 has two N-terminal entries; the isobaric filter keeps the one carrying diagnostic ions
    [TestCase(214, 'A', 144.102062, "iTRAQ-4plex on X")]
    public void TryResolve_NTerminalUnimodId_ResolvesToTheNTerminalMod(
        int unimodId, char? residue, double? mass, string expectedIdWithMotif)
    {
        var mod = CanonicalModification.AtNTerminus($"UNIMOD:{unimodId}", targetResidue: residue, mass: mass, unimodId: unimodId);

        var resolved = new UnimodModificationLookup().TryResolve(mod)?.MzLibModification;

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.IdWithMotif, Is.EqualTo(expectedIdWithMotif));
        Assert.That(resolved.LocationRestriction, Does.Contain("N-terminal"));
        Assert.That(UnimodIds(resolved), Does.Contain(unimodId.ToString()));
    }

    [TestCase(2, null, null, "Amidated on X")]
    [TestCase(2, 'K', null, "Amidated on X")]
    [TestCase(2, 'K', -0.984016, "Amidated on X")]
    [TestCase(34, 'K', null, "Methyl on X")]
    [TestCase(35, 'G', null, "Oxidation on G")]
    public void TryResolve_CTerminalUnimodId_ResolvesToTheCTerminalMod(
        int unimodId, char? residue, double? mass, string expectedIdWithMotif)
    {
        var mod = CanonicalModification.AtCTerminus($"UNIMOD:{unimodId}", targetResidue: residue, mass: mass, unimodId: unimodId);

        var resolved = new UnimodModificationLookup().TryResolve(mod)?.MzLibModification;

        Assert.That(resolved, Is.Not.Null);
        Assert.That(resolved!.IdWithMotif, Is.EqualTo(expectedIdWithMotif));
        Assert.That(resolved.LocationRestriction, Does.Contain("C-terminal"));
        Assert.That(UnimodIds(resolved), Does.Contain(unimodId.ToString()));
    }

    /// <summary>
    /// A terminal accession with no terminal specificity in Unimod still resolves, to the side-chain
    /// modification on that residue, as a terminal "Anywhere." modification always has.
    /// </summary>
    [TestCase(21, 'S', "Phospho on S")]
    [TestCase(35, 'M', "Oxidation on M")]
    public void TryResolve_NTerminalUnimodIdWithNoTerminalSpecificity_FallsBackToTheResidueMod(
        int unimodId, char residue, string expectedIdWithMotif)
    {
        var mod = CanonicalModification.AtNTerminus($"UNIMOD:{unimodId}", targetResidue: residue, unimodId: unimodId);

        var resolved = new UnimodModificationLookup().TryResolve(mod)?.MzLibModification;

        Assert.That(resolved?.IdWithMotif, Is.EqualTo(expectedIdWithMotif));
    }

    /// <summary>
    /// Amidation is C-terminal only in Unimod, so asking for it at the N-terminus has no answer.
    /// </summary>
    [Test]
    public void TryResolve_CTerminalOnlyUnimodIdAtTheNTerminus_IsUnresolved()
    {
        var mod = CanonicalModification.AtNTerminus("UNIMOD:2", targetResidue: 'K', unimodId: 2);

        var resolved = new UnimodModificationLookup().TryResolve(mod);

        Assert.That(resolved, Is.Null);
    }

    /// <summary>
    /// Residue modifications already resolved correctly; these pin that the terminal fix leaves them alone.
    /// </summary>
    [TestCase(21, 'S', 79.966331, "Phospho on S")]
    [TestCase(21, 'Y', null, "Phospho on Y")]
    [TestCase(35, 'M', 15.994915, "Oxidation on M")]
    [TestCase(4, 'C', null, "Carbamidomethyl on C")]
    [TestCase(1, 'K', null, "Acetyl on K")]
    [TestCase(121, 'K', null, "GG on K")]
    [TestCase(7, 'N', 0.984016, "Deamidated on N")]
    public void TryResolve_ResidueUnimodId_ResolvesToTheResidueMod(
        int unimodId, char residue, double? mass, string expectedIdWithMotif)
    {
        var mod = CanonicalModification.AtResidue(3, residue, $"UNIMOD:{unimodId}", mass: mass, unimodId: unimodId);

        var resolved = new UnimodModificationLookup().TryResolve(mod)?.MzLibModification;

        Assert.That(resolved?.IdWithMotif, Is.EqualTo(expectedIdWithMotif));
    }

    /// <summary>
    /// Met->Hse (UNIMOD:10) and Fluoro (UNIMOD:127) are on M and A and their ids contain "1", but
    /// Acetyl has no side-chain specificity on either residue.
    /// </summary>
    [TestCase('M')]
    [TestCase('A')]
    public void TryResolve_ResidueUnimodId_NeverResolvesToAModWithADifferentId(char residue)
    {
        var mod = CanonicalModification.AtResidue(3, residue, "UNIMOD:1", unimodId: 1);

        var resolved = new UnimodModificationLookup().TryResolve(mod)?.MzLibModification;

        Assert.That(resolved == null || UnimodIds(resolved).Contains("1"), Is.True,
            $"resolved to {resolved?.IdWithMotif} [{string.Join(",", UnimodIds(resolved!))}]");
    }

    #region Source modification's own Unimod reference (#1402)

    // A modification as ConvertModifications and ToCanonicalSequence hand it to the lookup: the mzLib id,
    // the formula and the Modification itself. Its mzLib id names nothing among the Unimod entries, so before
    // #1402 the formula fallback picked the shortest-named entry sharing the formula: Ethyl (UNIMOD:280) for
    // N6,N6-dimethyllysine, whose own reference is Dimethyl (UNIMOD:36).

    [TestCase("UniProt", "N6,N6-dimethyllysine on K", "Dimethyl on K", 36)]
    [TestCase("UniProt", "N6,N6,N6-trimethyllysine on K", "Trimethyl on K", 37)]
    [TestCase("UniProt", "Deamidated asparagine on N", "Deamidated on N", 7)]
    [TestCase("UniProt", "N-acetylserine on S", "Acetyl on S", 1)]
    [TestCase("Common Biological", "Dimethylation on K", "Dimethyl on K", 36)]
    [TestCase("Common Biological", "Trimethylation on K", "Trimethyl on K", 37)]
    [TestCase("Common Artifact", "Deamidation on N", "Deamidated on N", 7)]
    public void TryResolve_ModificationWithAUnimodReference_ResolvesToThatUnimodId(
        string modificationType, string idWithMotif, string expectedIdWithMotif, int expectedUnimodId)
    {
        var source = KnownMod(modificationType, idWithMotif);

        var resolved = new UnimodModificationLookup().TryResolve(AsConverted(source));

        Assert.That(resolved?.MzLibModification?.IdWithMotif, Is.EqualTo(expectedIdWithMotif));
        Assert.That(resolved!.Value.UnimodId, Is.EqualTo(expectedUnimodId));
        Assert.That(resolved.Value.MzLibModification!.ChemicalFormula, Is.EqualTo(source.ChemicalFormula));
    }

    /// <summary>
    /// The reference only chooses among entries with the source's formula. N,N-dimethylproline (C2H4) references
    /// Delta:H(5)C(2) (C2H5); taking it would change the mass, so the formula path keeps Dimethyl.
    /// </summary>
    [Test]
    public void TryResolve_UnimodReferenceWithADifferentFormula_KeepsTheFormulaMatch()
    {
        var source = KnownMod("UniProt", "N,N-dimethylproline on P");

        var resolved = new UnimodModificationLookup().TryResolve(AsConverted(source))?.MzLibModification;

        Assert.That(resolved?.IdWithMotif, Is.EqualTo("Dimethyl on P"));
        Assert.That(resolved!.ChemicalFormula, Is.EqualTo(source.ChemicalFormula));
    }

    /// <summary>
    /// With no Unimod reference the formula fallback and its shortest-name tie-break are unchanged: C2H4 on K
    /// is still Ethyl. Refusing an ambiguous formula match is out of scope for #1402.
    /// </summary>
    [Test]
    public void TryResolve_ModificationWithNoUnimodReference_KeepsTheFormulaFallback()
    {
        ModificationMotif.TryGetMotif("K", out var motifK);
        var source = new Modification(_originalId: "Custom C2H4", _modificationType: "Custom", _target: motifK,
            _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("C2H4"));

        var resolved = new UnimodModificationLookup().TryResolve(AsConverted(source))?.MzLibModification;

        Assert.That(resolved?.IdWithMotif, Is.EqualTo("Ethyl on K"));
    }

    /// <summary>
    /// MetaMorpheus modifications whose answer was already right, through the Unimod lookup, and every
    /// MetaMorpheus modification through the Global and mzLib lookups, which resolve it by mzLib id (Global may pick the same-named Unimod entry, as before).
    /// </summary>
    [TestCase("Common Variable", "Oxidation on M", "Oxidation on M")]
    [TestCase("Common Fixed", "Carbamidomethyl on C", "Carbamidomethyl on C")]
    [TestCase("Common Biological", "Phosphorylation on S", "Phospho on S")]
    [TestCase("Common Biological", "Acetylation on K", "Acetyl on K")]
    [TestCase("Common Biological", "Methylation on K", "Methyl on K")]
    public void TryResolve_MetaMorpheusModification_ResolvesAsBefore(
        string modificationType, string idWithMotif, string expectedUnimodIdWithMotif)
    {
        var source = KnownMod(modificationType, idWithMotif);

        Assert.That(new UnimodModificationLookup().TryResolve(AsConverted(source))?.MzLibModification?.IdWithMotif,
            Is.EqualTo(expectedUnimodIdWithMotif));
        Assert.That(new GlobalModificationLookup().TryResolve(AsConverted(source))?.MzLibModification?.IdWithMotif, Is.EqualTo(idWithMotif));
        Assert.That(new MzLibModificationLookup().TryResolve(AsConverted(source))?.MzLibModification, Is.SameAs(source));
    }

    /// <summary>
    /// The mzLib parser hands the lookup only the mzLib id, with no Modification and so no reference: unchanged.
    /// </summary>
    [Test]
    public void TryResolve_ParsedMzLibIdOnly_IsUnchanged()
    {
        const string representation = "UniProt:N6,N6-dimethyllysine on K";
        var mod = new CanonicalModification(ModificationPositionType.Residue, 3, 'K', representation, MzLibId: representation);

        Assert.That(new UnimodModificationLookup().TryResolve(mod), Is.Null);
    }

    [Test]
    public void ConvertModifications_UniProtDimethyllysine_BecomesUnimodDimethyl()
    {
        var dimethyl = KnownMod("UniProt", "N6,N6-dimethyllysine on K");
        var oxidation = KnownMod("Common Variable", "Oxidation on M");
        var protein = new Protein("PEPKMR", "TestProtein");
        var peptide = new PeptideWithSetModifications(protein, new DigestionParams(protease: "trypsin"),
            oneBasedStartResidueInProtein: 1, oneBasedEndResidueInProtein: 6, cleavageSpecificity: CleavageSpecificity.Full,
            peptideDescription: "Test", missedCleavages: 0,
            allModsOneIsNterminus: new Dictionary<int, Modification> { { 5, dimethyl }, { 6, oxidation } }, numFixedMods: 0);

        peptide.ConvertModifications(new UnimodModificationLookup());

        Assert.That(peptide.AllModsOneIsNterminus[5].IdWithMotif, Is.EqualTo("Dimethyl on K"));
        Assert.That(UnimodIds(peptide.AllModsOneIsNterminus[5]), Does.Contain("36"));
        Assert.That(peptide.AllModsOneIsNterminus[6].IdWithMotif, Is.EqualTo("Oxidation on M"));
    }

    /// <summary>
    /// A side-chain modification on a protein's first residue reaches the lookup as an N-terminal position. The
    /// Unimod reference must not move it to a terminal entry on another residue (Deamidated on F [N-terminal.]):
    /// the converted mod stays on its residue, so the tryptic digest still carries it.
    /// </summary>
    [TestCase("NPEPTIDEKAAAR", "Common Artifact", "Deamidation on N", "Deamidated on N")]
    [TestCase("RPEPTIDEKAAAR", "UniProt", "Citrulline on R", "Deamidated on R")]
    [TestCase("KPEPTIDEKAAAR", "Common Biological", "Carboxylation on K", "Carboxy on K")]
    public void ConvertModifications_SideChainModificationOnTheFirstResidue_StaysOnItsResidue(
        string sequence, string modificationType, string idWithMotif, string expectedIdWithMotif)
    {
        var source = KnownMod(modificationType, idWithMotif);
        var protein = new Protein(sequence, "TestProtein",
            oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { source } } });

        protein.ConvertModifications(new UnimodModificationLookup());

        var converted = protein.OneBasedPossibleLocalizedModifications[1].Single();
        Assert.That(converted.IdWithMotif, Is.EqualTo(expectedIdWithMotif));
        Assert.That(converted.Target.ToString(), Is.EqualTo(sequence[..1]));
        Assert.That(converted.ChemicalFormula, Is.EqualTo(source.ChemicalFormula));

        var firstPeptide = sequence[..(sequence.IndexOf('K', 1) + 1)];
        var modifiedPeptides = protein.Digest(new DigestionParams(protease: "trypsin"), new List<Modification>(), new List<Modification>())
            .Where(p => p.BaseSequence == firstPeptide && p.AllModsOneIsNterminus.Values.Contains(converted))
            .ToList();
        Assert.That(modifiedPeptides, Is.Not.Empty);
    }

    /// <summary>
    /// N-methylglycine on a mid-protein G references Methyl, but no Methyl entry restricted to a terminus may
    /// answer for a residue position (Methyl on X [Peptide C-terminal.]).
    /// </summary>
    [Test]
    public void TryResolve_NMethylglycineAtAResidue_IsNeverATerminalOnlyEntry()
    {
        var source = KnownMod("UniProt", "N-methylglycine on G");

        var resolved = new UnimodModificationLookup().TryResolve(AsConverted(source))?.MzLibModification;

        if (resolved != null)
        {
            Assert.That(resolved.LocationRestriction, Does.Not.Contain("terminal").IgnoreCase);
            Assert.That(resolved.Target.ToString(), Is.AnyOf("G", "X"));
            Assert.That(resolved.ChemicalFormula, Is.EqualTo(source.ChemicalFormula));
        }
    }

    private static Modification KnownMod(string modificationType, string idWithMotif) =>
        Mods.AllKnownMods.Single(m => m.ModificationType == modificationType && m.IdWithMotif == idWithMotif);

    /// <summary>What ConvertModifications builds for a residue modification.</summary>
    private static CanonicalModification AsConverted(Modification mod) =>
        CanonicalModification.AtResidue(3, mod.Target.ToString()[0], $"{mod.ModificationType}:{mod.IdWithMotif}",
            mod.MonoisotopicMass, mod.ChemicalFormula, mzLibId: mod.IdWithMotif, mzLibModification: mod);

    #endregion

    private static string[] UnimodIds(Modification modification) =>
        modification.DatabaseReference?
            .Where(kvp => kvp.Key.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase))
            .SelectMany(kvp => kvp.Value)
            .ToArray()
        ?? Array.Empty<string>();
}
