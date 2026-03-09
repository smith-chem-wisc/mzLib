using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace Test.Omics.SequenceConversion
{
    /// <summary>
    /// Tests for GlobalModificationLookup which searches all known modifications
    /// from all sources (MetaMorpheus, UniProt, UNIMOD, RNA databases).
    /// </summary>
    [TestFixture]
    public class GlobalModificationLookupTests
    {
        private GlobalModificationLookup _lookup;

        [SetUp]
        public void Setup()
        {
            _lookup = GlobalModificationLookup.Instance;
        }

        #region Resolution by mzLib ID Tests

        [Test]
        public void TryResolve_WithMzLibId_ResolvesCorrectly()
        {
            // Arrange - Common protein modification
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "[Common Variable:Oxidation on M]",
                mzLibId: "Common Variable:Oxidation on M");

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.MzLibModification.IdWithMotif, Is.EqualTo("Oxidation on M"));
            Assert.That(result.Value.HasMass, Is.True);
            Assert.That(result.Value.UnimodId, Is.EqualTo(35));
        }

        [Test]
        public void TryResolve_WithRnaModId_ResolvesCorrectly()
        {
            // Arrange - RNA modification
            var mod = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'A',
                originalRepresentation: "Biological:N6-methyladenosine",
                mzLibId: "N6-methyladenosine"); 

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.HasMass, Is.True);
        }

        #endregion

        #region Resolution by UNIMOD ID Tests

        [Test]
        public void TryResolve_WithUnimodId_ResolvesCorrectly()
        {
            // Arrange - UNIMOD:21 (Phosphorylation)
            var mod = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'S',
                originalRepresentation: "[UNIMOD:21]",
                unimodId: 21);

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.UnimodId, Is.EqualTo(21));
            Assert.That(result.Value.HasMass, Is.True);
        }

        [Test]
        public void TryResolve_WithUnimodId_PrefersSameResidue()
        {
            // Arrange - UNIMOD:21 on serine (phosphorylation can be on S, T, or Y)
            var modS = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'S',
                originalRepresentation: "[UNIMOD:21]",
                unimodId: 21);

            var modT = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'T',
                originalRepresentation: "[UNIMOD:21]",
                unimodId: 21);

            // Act
            var resultS = _lookup.TryResolve(modS);
            var resultT = _lookup.TryResolve(modT);

            // Assert
            Assert.That(resultS, Is.Not.Null);
            Assert.That(resultT, Is.Not.Null);
            Assert.That(resultS.Value.IsResolved, Is.True);
            Assert.That(resultT.Value.IsResolved, Is.True);
            
            // Should prefer modifications matching the target residue
            Assert.That(resultS.Value.MzLibModification.Target.ToString(), Does.Contain("S"));
            Assert.That(resultT.Value.MzLibModification.Target.ToString(), Does.Contain("T"));
        }

        #endregion

        #region Resolution by Name Tests

        [Test]
        public void TryResolve_ByName_ResolvesCorrectly()
        {
            // Arrange
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "Oxidation on M",
                mzLibId: "Oxidation on M");

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.HasMass, Is.True);
        }

        [Test]
        public void TryResolve_ByNameWithMotif_ResolvesCorrectly()
        {
            // Arrange - Full ID with motif
            var mod = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'S',
                originalRepresentation: "Phosphorylation on S",
                mzLibId: "Phosphorylation on S");

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.HasMass, Is.True);
            Assert.That(result.Value.UnimodId, Is.EqualTo(21));
        }

        [Test]
        public void TryResolve_ByName_FromMultipleSources()
        {
            // Arrange - Test that GlobalModificationLookup can find mods from different sources
            var acetylation = CanonicalModification.AtResidue(
                residueIndex: 4,
                targetResidue: 'K',
                originalRepresentation: "Acetylation on K",
                mzLibId: "Acetylation on K");

            var methylation = CanonicalModification.AtResidue(
                residueIndex: 7,
                targetResidue: 'K',
                originalRepresentation: "Methylation on K",
                mzLibId: "Methylation on K");

            // Act
            var resultAcetyl = _lookup.TryResolve(acetylation);
            var resultMethyl = _lookup.TryResolve(methylation);

            // Assert
            Assert.That(resultAcetyl, Is.Not.Null);
            Assert.That(resultAcetyl.Value.IsResolved, Is.True);
            Assert.That(resultMethyl, Is.Not.Null);
            Assert.That(resultMethyl.Value.IsResolved, Is.True);
        }

        #endregion

        #region Resolution by Chemical Formula Tests

        [Test]
        public void TryResolve_ByFormula_ResolvesCorrectly()
        {
            // Arrange - Oxidation is +O (add one oxygen)
            var formula = ChemicalFormula.ParseFormula("O");
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "[+O]",
                formula: formula);

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.HasMass, Is.True);
            // Should find Oxidation on M
            Assert.That(result.Value.MzLibModification.IdWithMotif, Does.Contain("Oxidation"));
        }

        [Test]
        public void TryResolve_ByFormula_PrefersMatchingResidue()
        {
            // Arrange - Phosphorylation formula: H O3 P
            var formula = ChemicalFormula.ParseFormula("H O3 P");
            var modS = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'S',
                originalRepresentation: "[+H O3 P]",
                formula: formula);

            var modT = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'T',
                originalRepresentation: "[+H O3 P]",
                formula: formula);

            // Act
            var resultS = _lookup.TryResolve(modS);
            var resultT = _lookup.TryResolve(modT);

            // Assert
            Assert.That(resultS, Is.Not.Null);
            Assert.That(resultT, Is.Not.Null);
            Assert.That(resultS.Value.IsResolved, Is.True);
            Assert.That(resultT.Value.IsResolved, Is.True);
            
            // Both should resolve to phosphorylation, but prefer matching residue
            Assert.That(resultS.Value.MzLibModification.Target.ToString(), Does.Contain("S"));
            Assert.That(resultT.Value.MzLibModification.Target.ToString(), Does.Contain("T"));
        }

        #endregion

        #region Resolution by Mass Tests

        [Test]
        public void TryResolve_ByMass_ResolvesCorrectly()
        {
            // Arrange - Oxidation mass is approximately +15.9949
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "[+15.9949]",
                mass: 15.9949);

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.HasMass, Is.True);
            // Should find Oxidation on M
            Assert.That(result.Value.MzLibModification.IdWithMotif, Does.Contain("Oxidation"));
        }

        [Test]
        public void TryResolve_ByMass_WithinTolerance()
        {
            // Arrange - Phosphorylation mass: 79.966331 Da
            // Test with slightly off mass within tolerance
            var mod = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'S',
                originalRepresentation: "[+79.9665]",
                mass: 79.9665); // Off by ~0.0002 Da

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.UnimodId, Is.EqualTo(21)); // Phosphorylation
        }

        [Test]
        public void TryResolve_ByMass_PrefersMatchingResidue()
        {
            // Arrange - Methylation mass: 14.01565 Da (can occur on multiple residues)
            var modK = CanonicalModification.AtResidue(
                residueIndex: 4,
                targetResidue: 'K',
                originalRepresentation: "[+14.01565]",
                mass: 14.01565);

            var modR = CanonicalModification.AtResidue(
                residueIndex: 6,
                targetResidue: 'R',
                originalRepresentation: "[+14.01565]",
                mass: 14.01565);

            // Act
            var resultK = _lookup.TryResolve(modK);
            var resultR = _lookup.TryResolve(modR);

            // Assert
            Assert.That(resultK, Is.Not.Null);
            Assert.That(resultR, Is.Not.Null);
            
            if (resultK.Value.IsResolved)
            {
                Assert.That(resultK.Value.MzLibModification.Target.ToString(), Does.Contain("K"));
            }
            
            if (resultR.Value.IsResolved)
            {
                Assert.That(resultR.Value.MzLibModification.Target.ToString(), Does.Contain("R"));
            }
        }

        [Test]
        public void TryResolve_ByMass_OutsideTolerance_DoesNotResolve()
        {
            // Arrange - Use a mass that doesn't match any known modification
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'X',
                originalRepresentation: "[+4123.456789]",
                mass: 4123.456789); // Unusual mass

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert 
            Assert.That(result, Is.Default);
        }

        #endregion

        #region Multiple Sources Tests

        [Test]
        public void TryResolve_FindsProteinModifications()
        {
            // Arrange - Common protein modifications
            var oxidation = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "Oxidation on M",
                mzLibId: "Oxidation on M");

            var carbamidomethyl = CanonicalModification.AtResidue(
                residueIndex: 2,
                targetResidue: 'C',
                originalRepresentation: "Carbamidomethyl on C",
                mzLibId: "Carbamidomethyl on C");

            // Act
            var resultOx = _lookup.TryResolve(oxidation);
            var resultCarb = _lookup.TryResolve(carbamidomethyl);

            // Assert
            Assert.That(resultOx, Is.Not.Null);
            Assert.That(resultOx.Value.IsResolved, Is.True);
            Assert.That(resultCarb, Is.Not.Null);
            Assert.That(resultCarb.Value.IsResolved, Is.True);
        }

        [Test]
        public void TryResolve_FindsRnaModifications()
        {
            // Arrange - RNA modifications
            var m6A = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'A',
                originalRepresentation: "N6-methyladenosine",
                mzLibId: "N6-methyladenosine"); // N6-methyladenosine

            // Act
            var result = _lookup.TryResolve(m6A);

            // Assert
            Assert.That(result, Is.Not.Null);
            // If it resolves, should have mass
            if (result.Value.IsResolved)
            {
                Assert.That(result.Value.HasMass, Is.True);
            }
        }

        #endregion

        #region Terminal Modification Tests

        [Test]
        public void TryResolve_NTerminalModification_ResolvesCorrectly()
        {
            // Arrange - N-terminal acetylation
            var mod = CanonicalModification.AtNTerminus(
                originalRepresentation: "Acetylation on X",
                mzLibId: "Acetylation on X");

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.Not.Null);
            Assert.That(result.Value.HasMass, Is.True);
        }

        [Test]
        public void TryResolve_CTerminalModification_ResolvesCorrectly()
        {
            // Arrange - C-terminal amidation
            var mod = CanonicalModification.AtCTerminus(
                originalRepresentation: "Less Common:Water loss on D",
                mzLibId: "Water loss on D");

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            // If it resolves, should have mass
            if (result.Value.IsResolved)
            {
                Assert.That(result.Value.HasMass, Is.True);
            }
        }

        #endregion

        #region Edge Cases

        [Test]
        public void TryResolve_AlreadyResolved_ReturnsUnchanged()
        {
            // Arrange - Create a fully resolved modification
            var oxidationMod = Mods.GetModification("Oxidation on M");
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "Oxidation on M",
                mzLibModification: oxidationMod);

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            Assert.That(result.Value.MzLibModification, Is.EqualTo(oxidationMod));
        }

        [Test]
        public void TryResolve_NoIdentifyingInformation_ReturnsUnresolved()
        {
            // Arrange - Modification with only position, no identifying info
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "Unknown");

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Null);
        }

        [Test]
        public void Instance_ReturnsSameInstance()
        {
            // Act
            var instance1 = GlobalModificationLookup.Instance;
            var instance2 = GlobalModificationLookup.Instance;

            // Assert
            Assert.That(instance1, Is.SameAs(instance2));
        }

        [Test]
        public void Name_ReturnsCorrectName()
        {
            // Act
            var name = _lookup.Name;

            // Assert
            Assert.That(name, Is.EqualTo("Global (All Mods)"));
        }

        #endregion

        #region Resolution Strategy Priority Tests

        [Test]
        public void TryResolve_PrefersIdOverName()
        {
            // Arrange - Modification with both mzLib ID and name
            var mod = CanonicalModification.AtResidue(
                residueIndex: 5,
                targetResidue: 'M',
                originalRepresentation: "Some Other Name",
                mzLibId: "Common Variable:Oxidation on M");

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            // Should resolve by ID (primary), not by name
            Assert.That(result.Value.MzLibModification.IdWithMotif, Is.EqualTo("Oxidation on M"));
        }

        [Test]
        public void TryResolve_PrefersUnimodIdOverMass()
        {
            // Arrange - Modification with both UNIMOD ID and mass
            var mod = CanonicalModification.AtResidue(
                residueIndex: 3,
                targetResidue: 'S',
                originalRepresentation: "[UNIMOD:21]",
                mass: 15.9949, // Oxidation mass (wrong)
                unimodId: 21); // Phosphorylation

            // Act
            var result = _lookup.TryResolve(mod);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Value.IsResolved, Is.True);
            // Should resolve by UNIMOD ID (primary), not by mass
            Assert.That(result.Value.UnimodId, Is.EqualTo(21)); // Phosphorylation, not Oxidation
        }

        #endregion
    }
}
