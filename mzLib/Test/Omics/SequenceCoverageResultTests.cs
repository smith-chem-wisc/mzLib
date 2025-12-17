using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SequenceCoverageResultTests
    {
        #region Constructor Tests

        [Test]
        public void Constructor_InitializesEmptyLists()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            Assert.That(result.SequenceCoverageFraction, Is.Not.Null);
            Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(0));

            Assert.That(result.SequenceCoverageDisplayList, Is.Not.Null);
            Assert.That(result.SequenceCoverageDisplayList.Count, Is.EqualTo(0));

            Assert.That(result.SequenceCoverageDisplayListWithMods, Is.Not.Null);
            Assert.That(result.SequenceCoverageDisplayListWithMods.Count, Is.EqualTo(0));

            Assert.That(result.FragmentSequenceCoverageDisplayList, Is.Not.Null);
            Assert.That(result.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(0));

            Assert.That(result.ModsInfo, Is.Not.Null);
            Assert.That(result.ModsInfo.Count, Is.EqualTo(0));
        }

        [Test]
        public void Constructor_ListsAreIndependent()
        {
            var result1 = new BioPolymerGroup.SequenceCoverageResult();
            var result2 = new BioPolymerGroup.SequenceCoverageResult();

            result1.SequenceCoverageFraction.Add(0.5);
            result1.SequenceCoverageDisplayList.Add("TEST");

            Assert.That(result2.SequenceCoverageFraction.Count, Is.EqualTo(0));
            Assert.That(result2.SequenceCoverageDisplayList.Count, Is.EqualTo(0));
        }

        #endregion

        #region SequenceCoverageFraction Tests

        [Test]
        public void SequenceCoverageFraction_CanAddValues()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.0);
            result.SequenceCoverageFraction.Add(0.5);
            result.SequenceCoverageFraction.Add(1.0);

            Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(3));
            Assert.That(result.SequenceCoverageFraction[0], Is.EqualTo(0.0));
            Assert.That(result.SequenceCoverageFraction[1], Is.EqualTo(0.5));
            Assert.That(result.SequenceCoverageFraction[2], Is.EqualTo(1.0));
        }

        [Test]
        public void SequenceCoverageFraction_AcceptsDecimalValues()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.12345);
            result.SequenceCoverageFraction.Add(0.99999);

            Assert.That(result.SequenceCoverageFraction[0], Is.EqualTo(0.12345).Within(0.00001));
            Assert.That(result.SequenceCoverageFraction[1], Is.EqualTo(0.99999).Within(0.00001));
        }

        [Test]
        public void SequenceCoverageFraction_CanClear()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.5);
            result.SequenceCoverageFraction.Add(0.75);
            result.SequenceCoverageFraction.Clear();

            Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(0));
        }

        [Test]
        public void SequenceCoverageFraction_CanContainDuplicates()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.5);
            result.SequenceCoverageFraction.Add(0.5);
            result.SequenceCoverageFraction.Add(0.5);

            Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(3));
        }

        #endregion

        #region SequenceCoverageDisplayList Tests

        [Test]
        public void SequenceCoverageDisplayList_CanAddStrings()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayList.Add("acdefghik");
            result.SequenceCoverageDisplayList.Add("ACDEFGHIK");

            Assert.That(result.SequenceCoverageDisplayList.Count, Is.EqualTo(2));
            Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo("acdefghik"));
            Assert.That(result.SequenceCoverageDisplayList[1], Is.EqualTo("ACDEFGHIK"));
        }

        [Test]
        public void SequenceCoverageDisplayList_AcceptsMixedCaseStrings()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayList.Add("acDEFGhik");

            Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo("acDEFGhik"));
        }

        [Test]
        public void SequenceCoverageDisplayList_AcceptsEmptyStrings()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayList.Add("");

            Assert.That(result.SequenceCoverageDisplayList.Count, Is.EqualTo(1));
            Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo(string.Empty));
        }

        [Test]
        public void SequenceCoverageDisplayList_AcceptsLongSequences()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();
            var longSequence = new string('A', 10000);

            result.SequenceCoverageDisplayList.Add(longSequence);

            Assert.That(result.SequenceCoverageDisplayList[0].Length, Is.EqualTo(10000));
        }

        #endregion

        #region SequenceCoverageDisplayListWithMods Tests

        [Test]
        public void SequenceCoverageDisplayListWithMods_CanAddModifiedSequences()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayListWithMods.Add("[Acetyl on A]-ACDEFghik");
            result.SequenceCoverageDisplayListWithMods.Add("acD[Phospho on S]EFghik");

            Assert.That(result.SequenceCoverageDisplayListWithMods.Count, Is.EqualTo(2));
            Assert.That(result.SequenceCoverageDisplayListWithMods[0], Does.Contain("[Acetyl on A]"));
            Assert.That(result.SequenceCoverageDisplayListWithMods[1], Does.Contain("[Phospho on S]"));
        }

        [Test]
        public void SequenceCoverageDisplayListWithMods_AcceptsNTerminalMods()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayListWithMods.Add("[Acetyl on A]-ACDEFGHIK");

            Assert.That(result.SequenceCoverageDisplayListWithMods[0], Does.StartWith("[Acetyl on A]-"));
        }

        [Test]
        public void SequenceCoverageDisplayListWithMods_AcceptsCTerminalMods()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayListWithMods.Add("ACDEFGHIK-[Amidation on K]");

            Assert.That(result.SequenceCoverageDisplayListWithMods[0], Does.EndWith("-[Amidation on K]"));
        }

        [Test]
        public void SequenceCoverageDisplayListWithMods_AcceptsMultipleMods()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayListWithMods.Add("[Acetyl on A]-ACD[Phospho on S]EFghik-[Amidation on K]");

            var entry = result.SequenceCoverageDisplayListWithMods[0];
            Assert.That(entry, Does.Contain("[Acetyl on A]"));
            Assert.That(entry, Does.Contain("[Phospho on S]"));
            Assert.That(entry, Does.Contain("[Amidation on K]"));
        }

        #endregion

        #region FragmentSequenceCoverageDisplayList Tests

        [Test]
        public void FragmentSequenceCoverageDisplayList_CanAddStrings()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.FragmentSequenceCoverageDisplayList.Add("acdefghik");
            result.FragmentSequenceCoverageDisplayList.Add("ACDefghik");

            Assert.That(result.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(2));
        }

        [Test]
        public void FragmentSequenceCoverageDisplayList_IndependentFromSequenceCoverageDisplayList()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayList.Add("ACDEFGHIK");
            result.FragmentSequenceCoverageDisplayList.Add("acdefghik");

            Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo("ACDEFGHIK"));
            Assert.That(result.FragmentSequenceCoverageDisplayList[0], Is.EqualTo("acdefghik"));
            Assert.That(result.SequenceCoverageDisplayList.Count, Is.EqualTo(1));
            Assert.That(result.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(1));
        }

        [Test]
        public void FragmentSequenceCoverageDisplayList_CanRepresentPartialCoverage()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            // Fragment coverage typically shows partial coverage
            result.FragmentSequenceCoverageDisplayList.Add("ACDefghik"); // Only first 3 residues covered

            var coverage = result.FragmentSequenceCoverageDisplayList[0];
            Assert.That(coverage.Substring(0, 3), Is.EqualTo("ACD")); // Covered (uppercase)
            Assert.That(coverage.Substring(3), Is.EqualTo("efghik")); // Not covered (lowercase)
        }

        #endregion

        #region ModsInfo Tests

        [Test]
        public void ModsInfo_CanAddOccupancyStrings()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.ModsInfo.Add("#aa3[Phospho on S,info:occupancy=0.50(1/2)]");

            Assert.That(result.ModsInfo.Count, Is.EqualTo(1));
            Assert.That(result.ModsInfo[0], Does.Contain("#aa3"));
            Assert.That(result.ModsInfo[0], Does.Contain("occupancy"));
        }

        [Test]
        public void ModsInfo_CanAddMultipleModifications()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.ModsInfo.Add("#aa3[Phospho on S,info:occupancy=0.50(1/2)];#aa7[Acetyl on K,info:occupancy=1.00(2/2)]");

            Assert.That(result.ModsInfo[0], Does.Contain("#aa3"));
            Assert.That(result.ModsInfo[0], Does.Contain("#aa7"));
        }

        [Test]
        public void ModsInfo_CanAddMultipleEntriesForDifferentProteins()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.ModsInfo.Add("#aa5[Phospho on S,info:occupancy=1.00(3/3)]");
            result.ModsInfo.Add("#aa10[Oxidation on M,info:occupancy=0.33(1/3)]");

            Assert.That(result.ModsInfo.Count, Is.EqualTo(2));
        }

        [Test]
        public void ModsInfo_AcceptsEmptyString()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.ModsInfo.Add("");

            Assert.That(result.ModsInfo.Count, Is.EqualTo(1));
            Assert.That(result.ModsInfo[0], Is.EqualTo(string.Empty));
        }

        [Test]
        public void ModsInfo_OccupancyFormatIsCorrect()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            // Expected format: #aa{position}[{modName},info:occupancy={fraction}({count}/{total})]
            var modInfo = "#aa15[Phosphorylation on S,info:occupancy=0.75(3/4)]";
            result.ModsInfo.Add(modInfo);

            Assert.That(result.ModsInfo[0], Does.Match(@"#aa\d+\[.+,info:occupancy=\d+\.\d+\(\d+/\d+\)\]"));
        }

        #endregion

        #region List Behavior Tests

        [Test]
        public void AllLists_SupportAddRange()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.AddRange(new[] { 0.1, 0.2, 0.3 });
            result.SequenceCoverageDisplayList.AddRange(new[] { "SEQ1", "SEQ2" });
            result.SequenceCoverageDisplayListWithMods.AddRange(new[] { "MOD1", "MOD2" });
            result.FragmentSequenceCoverageDisplayList.AddRange(new[] { "FRAG1", "FRAG2" });
            result.ModsInfo.AddRange(new[] { "INFO1", "INFO2" });

            Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(3));
            Assert.That(result.SequenceCoverageDisplayList.Count, Is.EqualTo(2));
            Assert.That(result.SequenceCoverageDisplayListWithMods.Count, Is.EqualTo(2));
            Assert.That(result.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(2));
            Assert.That(result.ModsInfo.Count, Is.EqualTo(2));
        }

        [Test]
        public void AllLists_SupportRemove()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.5);
            result.SequenceCoverageFraction.Remove(0.5);

            result.SequenceCoverageDisplayList.Add("TEST");
            result.SequenceCoverageDisplayList.Remove("TEST");

            Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(0));
            Assert.That(result.SequenceCoverageDisplayList.Count, Is.EqualTo(0));
        }

        [Test]
        public void AllLists_SupportIndexAccess()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.25);
            result.SequenceCoverageFraction.Add(0.75);

            Assert.That(result.SequenceCoverageFraction[0], Is.EqualTo(0.25));
            Assert.That(result.SequenceCoverageFraction[1], Is.EqualTo(0.75));
        }

        [Test]
        public void AllLists_SupportContains()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageDisplayList.Add("ACDEFGHIK");

            Assert.That(result.SequenceCoverageDisplayList.Contains("ACDEFGHIK"), Is.True);
            Assert.That(result.SequenceCoverageDisplayList.Contains("NOTFOUND"), Is.False);
        }

        #endregion

        #region Integration with BioPolymerGroup Tests

        [Test]
        public void SequenceCoverageResult_UsedByBioPolymerGroupToString()
        {
            // Create a minimal BioPolymerGroup setup
            var bioPolymer = new TestBioPolymer("ACDEFGHIK", "P00001");
            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods>();
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods>();

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);

            // Before CalculateSequenceCoverage, ToString should not throw
            Assert.DoesNotThrow(() => group.ToString());

            // After CalculateSequenceCoverage, ToString should include coverage data
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>();
            group.CalculateSequenceCoverage();

            var output = group.ToString();
            Assert.That(output, Is.Not.Null);
            Assert.That(output, Is.Not.Empty);
        }

        [Test]
        public void SequenceCoverageResult_PopulatedByCalculateSequenceCoverage()
        {
            var bioPolymer = new TestBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", parent: bioPolymer, startResidue: 1, endResidue: 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            // Verify coverage is reflected in ToString output
            var output = group.ToString();

            // Coverage fraction should be 5/9 ≈ 0.556
            Assert.That(output, Does.Contain("0.5")); // Partial match on fraction

            // Coverage display should show uppercase for covered residues
            Assert.That(output, Does.Contain("ACDEF")); // First 5 residues covered
        }

        [Test]
        public void SequenceCoverageResult_ReplacedOnMultipleCalculations()
        {
            var bioPolymer = new TestBioPolymer("ACDEFGHIK", "P00001");
            var peptide = new CoverageBioPolymerWithSetMods("ACDEF", "ACDEF", parent: bioPolymer, startResidue: 1, endResidue: 5);

            var bioPolymers = new HashSet<IBioPolymer> { bioPolymer };
            var sequences = new HashSet<IBioPolymerWithSetMods> { peptide };
            var uniqueSequences = new HashSet<IBioPolymerWithSetMods> { peptide };

            var psm = new BioPolymerGroupSequenceCoverageTests.CoverageSpectralMatch(@"C:\test.raw", "ACDEF", "ACDEF", 100, 1, new[] { peptide });

            var group = new BioPolymerGroup(bioPolymers, sequences, uniqueSequences);
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();
            var output1 = group.ToString();

            group.CalculateSequenceCoverage();
            var output2 = group.ToString();

            // Results should be identical (replaced, not appended)
            Assert.That(output1, Is.EqualTo(output2));
        }

        #endregion

        #region Edge Cases

        [Test]
        public void SequenceCoverageResult_HandlesSpecialCharactersInModNames()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.ModsInfo.Add("#aa5[Phospho (STY),info:occupancy=0.50(1/2)]");
            result.SequenceCoverageDisplayListWithMods.Add("acde[Phospho (STY)]fghik");

            Assert.That(result.ModsInfo[0], Does.Contain("Phospho (STY)"));
            Assert.That(result.SequenceCoverageDisplayListWithMods[0], Does.Contain("[Phospho (STY)]"));
        }

        [Test]
        public void SequenceCoverageResult_HandlesSingleResidueProtein()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(1.0);
            result.SequenceCoverageDisplayList.Add("K");
            result.FragmentSequenceCoverageDisplayList.Add("K");

            Assert.That(result.SequenceCoverageFraction[0], Is.EqualTo(1.0));
            Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo("K"));
        }

        [Test]
        public void SequenceCoverageResult_HandlesVeryLongProtein()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            var longSequence = new string('a', 10000);
            result.SequenceCoverageDisplayList.Add(longSequence);
            result.SequenceCoverageFraction.Add(0.001);

            Assert.That(result.SequenceCoverageDisplayList[0].Length, Is.EqualTo(10000));
            Assert.That(result.SequenceCoverageFraction[0], Is.EqualTo(0.001));
        }

        [Test]
        public void SequenceCoverageResult_HandlesZeroCoverage()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.0);
            result.SequenceCoverageDisplayList.Add("acdefghik"); // All lowercase = no coverage

            Assert.That(result.SequenceCoverageFraction[0], Is.EqualTo(0.0));
            Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo(result.SequenceCoverageDisplayList[0].ToLower()));
        }

        [Test]
        public void SequenceCoverageResult_HandlesFullCoverage()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(1.0);
            result.SequenceCoverageDisplayList.Add("ACDEFGHIK"); // All uppercase = full coverage

            Assert.That(result.SequenceCoverageFraction[0], Is.EqualTo(1.0));
            Assert.That(result.SequenceCoverageDisplayList[0], Is.EqualTo(result.SequenceCoverageDisplayList[0].ToUpper()));
        }

        [Test]
        public void SequenceCoverageResult_ListsCanBeEnumerated()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            result.SequenceCoverageFraction.Add(0.1);
            result.SequenceCoverageFraction.Add(0.2);
            result.SequenceCoverageFraction.Add(0.3);

            var sum = 0.0;
            foreach (var fraction in result.SequenceCoverageFraction)
            {
                sum += fraction;
            }

            Assert.That(sum, Is.EqualTo(0.6).Within(0.001));
        }

        [Test]
        public void SequenceCoverageResult_MultipleProteinsInGroup()
        {
            var result = new BioPolymerGroup.SequenceCoverageResult();

            // Simulate a group with 3 proteins
            result.SequenceCoverageFraction.Add(0.3);
            result.SequenceCoverageFraction.Add(0.5);
            result.SequenceCoverageFraction.Add(0.8);

            result.SequenceCoverageDisplayList.Add("ACDEfghik");
            result.SequenceCoverageDisplayList.Add("MNPQRstvwy");
            result.SequenceCoverageDisplayList.Add("ABCDEFGH");

            result.FragmentSequenceCoverageDisplayList.Add("acdefghik");
            result.FragmentSequenceCoverageDisplayList.Add("mnpqrstvwy");
            result.FragmentSequenceCoverageDisplayList.Add("abcdefgh");

            Assert.That(result.SequenceCoverageFraction.Count, Is.EqualTo(3));
            Assert.That(result.SequenceCoverageDisplayList.Count, Is.EqualTo(3));
            Assert.That(result.FragmentSequenceCoverageDisplayList.Count, Is.EqualTo(3));
        }

        #endregion
    }

    #region Additional Test Helper Classes

    /// <summary>
    /// Extended test implementation of IBioPolymerWithSetMods that supports Parent property.
    /// Used specifically for SequenceCoverageResult tests.
    /// </summary>
    internal class CoverageBioPolymerWithSetMods : IBioPolymerWithSetMods
    {
        public string BaseSequence { get; }
        public string FullSequence { get; }
        public double MostAbundantMonoisotopicMass { get; }
        public string SequenceWithChemicalFormulas { get; }
        public int OneBasedStartResidue { get; }
        public int OneBasedEndResidue { get; }
        public int MissedCleavages { get; }
        public string Description { get; }
        public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; }
        public char PreviousResidue { get; }
        public char NextResidue { get; }
        public IDigestionParams DigestionParams { get; }
        public Dictionary<int, Modification> AllModsOneIsNterminus { get; }
        public int NumMods { get; }
        public int NumFixedMods { get; }
        public int NumVariableMods { get; }
        public int Length => BaseSequence.Length;
        public IBioPolymer Parent { get; }
        public ChemicalFormula ThisChemicalFormula { get; }
        public double MonoisotopicMass { get; }

        public CoverageBioPolymerWithSetMods(string baseSequence, string fullSequence,
            double mass = 0, int startResidue = 1, int endResidue = 0, IBioPolymer parent = null)
        {
            BaseSequence = baseSequence;
            FullSequence = fullSequence;
            MostAbundantMonoisotopicMass = mass;
            MonoisotopicMass = mass;
            SequenceWithChemicalFormulas = baseSequence;
            OneBasedStartResidue = startResidue;
            OneBasedEndResidue = endResidue > 0 ? endResidue : startResidue + baseSequence.Length - 1;
            MissedCleavages = 0;
            Description = "Test";
            CleavageSpecificityForFdrCategory = CleavageSpecificity.Full;
            PreviousResidue = '-';
            NextResidue = '-';
            DigestionParams = null;
            AllModsOneIsNterminus = new Dictionary<int, Modification>();
            NumMods = 0;
            NumFixedMods = 0;
            NumVariableMods = 0;
            Parent = parent;
            ThisChemicalFormula = new ChemicalFormula();
        }

        public void Fragment(DissociationType dissociationType, FragmentationTerminus fragmentationTerminus,
            List<Product> products, FragmentationParams? fragmentationParams = null)
        {
        }

        public void FragmentInternally(DissociationType dissociationType, int minLengthOfFragments,
            List<Product> products, FragmentationParams? fragmentationParams = null)
        {
        }

        public IBioPolymerWithSetMods Localize(int indexOfMass, double massToLocalize)
        {
            return this;
        }

        public bool Equals(IBioPolymerWithSetMods? other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            return BaseSequence == other.BaseSequence && FullSequence == other.FullSequence;
        }

        public override bool Equals(object? obj)
        {
            return obj is IBioPolymerWithSetMods other && Equals(other);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(BaseSequence, FullSequence);
        }
    }

    #endregion
}