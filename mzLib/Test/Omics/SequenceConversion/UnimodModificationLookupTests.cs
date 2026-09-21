using System;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

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

    private static string[] UnimodIds(Modification modification) =>
        modification.DatabaseReference?
            .Where(kvp => kvp.Key.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase))
            .SelectMany(kvp => kvp.Value)
            .ToArray()
        ?? Array.Empty<string>();
}
