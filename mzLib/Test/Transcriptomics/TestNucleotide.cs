using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Chemistry;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace Test.Transcriptomics
{
    [ExcludeFromCodeCoverage]
    [TestFixture]
    public class TestNucleotide
    {
        public record NucleotideTestCase(Nucleotide Nucleotide, string Name, char OneLetterCode, string Symbol, ChemicalFormula Formula, double Mass,
            ChemicalFormula nucleosideFormula);

        public static IEnumerable<NucleotideTestCase> GetNucleotideTestCases()
        {
            yield return new NucleotideTestCase(Nucleotide.AdenineBase, "Adenine", 'A', "Ade",
                ChemicalFormula.ParseFormula("C5H4N5"), 329.052523, ChemicalFormula.ParseFormula("C10H13N5O4"));
            yield return new NucleotideTestCase(Nucleotide.CytosineBase, "Cytosine", 'C',
                "Cyt", ChemicalFormula.ParseFormula("C4H4N3O1"), 305.041290, ChemicalFormula.ParseFormula("C9H13N3O5"));
            yield return new NucleotideTestCase(Nucleotide.GuanineBase, "Guanine", 'G',
                "Gua", ChemicalFormula.ParseFormula("C5H4N5O1"), 345.047438, ChemicalFormula.ParseFormula("C10H13N5O5"));
            yield return new NucleotideTestCase(Nucleotide.UracilBase, "Uracil", 'U',
                "Ura", ChemicalFormula.ParseFormula("C4H3N2O2"), 306.025306, ChemicalFormula.ParseFormula("C9H12N2O6"));
            yield return new NucleotideTestCase(Nucleotide.DeoxyAdenineBase, "DeoxyAdenine", 'B',
                "dAde", ChemicalFormula.ParseFormula("C5H4N5"), 313.057607, ChemicalFormula.ParseFormula("C10H13N5O3"));
            yield return new NucleotideTestCase(Nucleotide.DeoxyCytosineBase, "DeoxyCytosine", 'D',
                "dCyt", ChemicalFormula.ParseFormula("C4H4N3O1"), 289.046375, ChemicalFormula.ParseFormula("C9H13N3O4"));
            yield return new NucleotideTestCase(Nucleotide.DeoxyGuanineBase, "DeoxyGuanine", 'H',
                "dGua", ChemicalFormula.ParseFormula("C5H4N5O1"), 329.052523, ChemicalFormula.ParseFormula("C10H13N5O4"));
            yield return new NucleotideTestCase(Nucleotide.DeoxyThymineBase, "DeoxyThymine", 'V',
                "dThy", ChemicalFormula.ParseFormula("C5H5N2O2"), 304.046041, ChemicalFormula.ParseFormula("C10H14N2O5"));
            yield return new NucleotideTestCase(Nucleotide.PseudoUracilBase, "PseudoUracil", 'Y',
                "Psu", ChemicalFormula.ParseFormula("C4H3N2O2"), 306.025306, ChemicalFormula.ParseFormula("C9H12N2O6"));
        }

        [Test]
        [TestCaseSource(nameof(GetNucleotideTestCases))]
        public void TestCommonNucleotides(NucleotideTestCase testCase)
        {
            Nucleotide nucleotide = testCase.Nucleotide;

            Assert.That(nucleotide.MonoisotopicMass, Is.EqualTo(testCase.Mass).Within(0.00001));
            Assert.That(nucleotide.Letter, Is.EqualTo(testCase.OneLetterCode));
            Assert.That(nucleotide.Symbol, Is.EqualTo(testCase.Symbol));
            Assert.That(nucleotide.ToString(), Is.EqualTo($"{testCase.OneLetterCode} {testCase.Symbol} ({testCase.Name})"));
            Assert.That(nucleotide.BaseChemicalFormula, Is.EqualTo(testCase.Formula));
            Assert.That(nucleotide.NucleosideChemicalFormula, Is.EqualTo(testCase.nucleosideFormula));

            Nucleotide newNucleotide =
                new Nucleotide(testCase.Name, testCase.OneLetterCode, testCase.Symbol, testCase.Formula);
            Assert.That(nucleotide.Equals(nucleotide));
            Assert.That(!nucleotide.Equals(null));
            Assert.That(nucleotide.Equals(newNucleotide));
            Assert.That(nucleotide.Equals((object)newNucleotide));
            Assert.That(!nucleotide.Equals((object)null));
        }

        [Test]
        [TestCaseSource(nameof(GetNucleotideTestCases))]
        public void TestGetResidue(NucleotideTestCase testCase)
        {
            Nucleotide nucleotide = testCase.Nucleotide;

            var testNucleotide = Nucleotide.GetResidue(testCase.OneLetterCode);
            Assert.That(nucleotide.Equals(testNucleotide));

            if (Nucleotide.TryGetResidue(testCase.OneLetterCode, out Nucleotide outTide))
            {
                Assert.That(nucleotide.Equals(outTide));
            }
            else
                Assert.Fail();

            testNucleotide = Nucleotide.GetResidue(testCase.Symbol);
            Assert.That(nucleotide.Equals(testNucleotide));
            if (Nucleotide.TryGetResidue(testCase.Symbol, out outTide))
            {
                Assert.That(nucleotide.Equals(outTide));
            }
            else
                Assert.Fail();

            if (Nucleotide.TryGetResidue('&', out outTide))
                Assert.Fail();
            else
                Assert.Pass();
        }

        [Test]
        public void TestPseudoUracilAlternatePsiLookup()
        {
            const char psi = '\u03A8';
            const char lowerPsi = '\u03C8';

            Assert.That(Nucleotide.TryGetResidue(psi, out Nucleotide psiTide), Is.True);
            Assert.That(psiTide, Is.EqualTo(Nucleotide.PseudoUracilBase));

            Assert.That(Nucleotide.GetResidue(psi), Is.EqualTo(Nucleotide.PseudoUracilBase));
            Assert.That(Nucleotide.GetResidue(psi.ToString()), Is.EqualTo(Nucleotide.PseudoUracilBase));

            Assert.That(Nucleotide.TryGetResidue(lowerPsi, out Nucleotide outTide), Is.False);
            Assert.That(outTide, Is.Null);
            Assert.That(Nucleotide.TryGetResidue(lowerPsi.ToString(), out outTide), Is.False);
            Assert.That(outTide, Is.Null);
            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue(lowerPsi));
            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue(lowerPsi.ToString()));

            Assert.That(Nucleotide.PseudoUracilBase.Letter, Is.EqualTo('Y'));
        }

        [Test]
        public void TestBuiltInResiduesAreExactCaseOnly()
        {
            Assert.That(Nucleotide.GetResidue('A'), Is.EqualTo(Nucleotide.AdenineBase));
            Assert.That(Nucleotide.GetResidue("A"), Is.EqualTo(Nucleotide.AdenineBase));

            Assert.That(Nucleotide.TryGetResidue('a', out Nucleotide outTide), Is.False);
            Assert.That(outTide, Is.Null);
            Assert.That(Nucleotide.TryGetResidue("a", out outTide), Is.False);
            Assert.That(outTide, Is.Null);

            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue('a'));
            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue("a"));
        }

        [Test]
        public void TestCustomResidueWithAlternateCodes()
        {
            string name = "FakeNucleotidePsi";
            char oneLetter = 'Q';
            string symbol = "Fkp";
            string chemicalFormula = "C5H5N2O2";
            char alternateChar = '\u03A9';
            string alternateString = "FkpAlt";
            var fakeNucleotide = new Nucleotide(name, oneLetter, symbol, ChemicalFormula.ParseFormula(chemicalFormula));

            Nucleotide.AddResidue(name, oneLetter, symbol, chemicalFormula);

            Assert.That(Nucleotide.TryAddAlternativeRepresentation(fakeNucleotide, alternateChar), Is.True);
            Assert.That(Nucleotide.TryAddAlternativeRepresentation(fakeNucleotide, alternateString), Is.True);

            Assert.That(Nucleotide.GetResidue(alternateChar), Is.EqualTo(fakeNucleotide));
            Assert.That(Nucleotide.TryGetResidue(alternateChar, out Nucleotide outTide), Is.True);
            Assert.That(outTide, Is.EqualTo(fakeNucleotide));
            Assert.That(Nucleotide.GetResidue(alternateString), Is.EqualTo(fakeNucleotide));

            Assert.That(Nucleotide.TryAddAlternativeRepresentation(fakeNucleotide, alternateChar), Is.True);
            Assert.That(Nucleotide.TryAddAlternativeRepresentation(fakeNucleotide, alternateString), Is.True);

            Assert.That(Nucleotide.GetResidue(oneLetter), Is.EqualTo(fakeNucleotide));
            Assert.That(fakeNucleotide.Letter, Is.EqualTo(oneLetter));
        }

        [Test]
        public void TestTryAddAlternativeRepresentationClash()
        {
            string name = "ClashNucleotide";
            char oneLetter = 'R';
            string symbol = "Cls";
            string chemicalFormula = "C5H5N2O3";
            var clashNucleotide = new Nucleotide(name, oneLetter, symbol, ChemicalFormula.ParseFormula(chemicalFormula));
            Nucleotide.AddResidue(name, oneLetter, symbol, chemicalFormula);

            Assert.That(Nucleotide.TryAddAlternativeRepresentation(clashNucleotide, "Psu"), Is.False);
            Assert.That(Nucleotide.GetResidue("Psu"), Is.EqualTo(Nucleotide.PseudoUracilBase));

            Assert.That(Nucleotide.TryAddAlternativeRepresentation(clashNucleotide, 'Y'), Is.False);
            Assert.That(Nucleotide.GetResidue('Y'), Is.EqualTo(Nucleotide.PseudoUracilBase));
        }

        [Test]
        public void TestTryAddAlternativeRepresentationAsciiCharClashLeavesLookupUnchanged()
        {
            string originalName = "ExistingAsciiAliasNucleotide";
            char originalOneLetter = 'L';
            string originalSymbol = "Eaa";
            string originalChemicalFormula = "C5H5N2O5";
            var existingNucleotide = new Nucleotide(originalName, originalOneLetter, originalSymbol,
                ChemicalFormula.ParseFormula(originalChemicalFormula));
            Nucleotide.AddResidue(originalName, originalOneLetter, originalSymbol, originalChemicalFormula);
            Assert.That(Nucleotide.TryAddAlternativeRepresentation(existingNucleotide, "J"), Is.True);

            string conflictingName = "ConflictingAsciiAliasNucleotide";
            char conflictingOneLetter = 'M';
            string conflictingSymbol = "Caa";
            string conflictingChemicalFormula = "C5H5N2O6";
            var conflictingNucleotide = new Nucleotide(conflictingName, conflictingOneLetter, conflictingSymbol,
                ChemicalFormula.ParseFormula(conflictingChemicalFormula));
            Nucleotide.AddResidue(conflictingName, conflictingOneLetter, conflictingSymbol, conflictingChemicalFormula);

            Assert.That(Nucleotide.TryAddAlternativeRepresentation(conflictingNucleotide, 'J'), Is.False);
            Assert.That(Nucleotide.GetResidue('J'), Is.EqualTo(existingNucleotide));
            Assert.That(Nucleotide.GetResidue("J"), Is.EqualTo(existingNucleotide));
            Assert.That(Nucleotide.TryGetResidue('J', out Nucleotide outTide), Is.True);
            Assert.That(outTide, Is.EqualTo(existingNucleotide));
            Assert.That(Nucleotide.TryGetResidue("J", out outTide), Is.True);
            Assert.That(outTide, Is.EqualTo(existingNucleotide));
        }

        [Test]
        public void TestPseudoUracilPsiSequenceRoundTrip()
        {
            const char psi = '\u03A8';
            var rna = new RNA($"GU{psi}CUG", "");

            Assert.That(rna.BaseSequence, Is.EqualTo("GUYCUG"));
        }

        [Test]
        public void TestNonAsciiUnknownCharLookup()
        {
            const char alpha = '\u03B1';

            Assert.That(Nucleotide.TryGetResidue(alpha, out Nucleotide outTide), Is.False);
            Assert.That(outTide, Is.Null);

            Assert.That(Nucleotide.TryGetResidue('&', out outTide), Is.False);
            Assert.That(Nucleotide.TryGetResidue("Taco", out outTide), Is.False);

            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue(alpha));
            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue('&'));
            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue("Taco"));
        }

        [Test]
        public void TestTryAddAlternativeRepresentationAsciiChar()
        {
            string name = "AsciiAliasNucleotide";
            char oneLetter = 'X';
            string symbol = "Asl";
            string chemicalFormula = "C5H5N2O4";
            var asciiNucleotide = new Nucleotide(name, oneLetter, symbol, ChemicalFormula.ParseFormula(chemicalFormula));
            Nucleotide.AddResidue(name, oneLetter, symbol, chemicalFormula);

            char alias = 'W';
            Assert.That(Nucleotide.TryAddAlternativeRepresentation(asciiNucleotide, alias), Is.True);
            Assert.That(Nucleotide.GetResidue(alias), Is.EqualTo(asciiNucleotide));
            Assert.That(Nucleotide.GetResidue(alias.ToString()), Is.EqualTo(asciiNucleotide));
            Assert.That(Nucleotide.TryGetResidue(alias, out Nucleotide outTide), Is.True);
            Assert.That(outTide, Is.EqualTo(asciiNucleotide));
            Assert.That(Nucleotide.TryGetResidue(alias.ToString(), out outTide), Is.True);
            Assert.That(outTide, Is.EqualTo(asciiNucleotide));

            Assert.That(Nucleotide.TryGetResidue('w', out outTide), Is.False);
            Assert.That(outTide, Is.Null);
            Assert.That(Nucleotide.TryGetResidue("w", out outTide), Is.False);
            Assert.That(outTide, Is.Null);
            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue('w'));
            Assert.Throws<KeyNotFoundException>(() => Nucleotide.GetResidue("w"));

            Assert.That(Nucleotide.TryAddAlternativeRepresentation(asciiNucleotide, alias), Is.True);
        }

        [Test]
        public void TestTryAddAlternativeRepresentationSingleCharacterStringUsesExactCaseLookup()
        {
            string name = "SingleCharacterStringAliasNucleotide";
            char oneLetter = 'N';
            string symbol = "Sca";
            string chemicalFormula = "C5H5N3O4";
            var stringAliasNucleotide = new Nucleotide(name, oneLetter, symbol, ChemicalFormula.ParseFormula(chemicalFormula));
            Nucleotide.AddResidue(name, oneLetter, symbol, chemicalFormula);

            Assert.That(Nucleotide.TryAddAlternativeRepresentation(stringAliasNucleotide, "Z"), Is.True);
            Assert.That(Nucleotide.TryGetResidue("Z", out Nucleotide outTide), Is.True);
            Assert.That(outTide, Is.EqualTo(stringAliasNucleotide));
            Assert.That(Nucleotide.TryGetResidue('Z', out outTide), Is.True);
            Assert.That(outTide, Is.EqualTo(stringAliasNucleotide));
            Assert.That(Nucleotide.GetResidue("Z"), Is.EqualTo(stringAliasNucleotide));
            Assert.That(Nucleotide.GetResidue('Z'), Is.EqualTo(stringAliasNucleotide));

            var rna = new RNA("GUZCUG", "");
            Assert.That(rna.BaseSequence, Is.EqualTo("GUNCUG"));

            Assert.That(Nucleotide.TryGetResidue("k", out outTide), Is.False);
            Assert.That(outTide, Is.Null);
            Assert.That(Nucleotide.TryGetResidue('k', out outTide), Is.False);
            Assert.That(outTide, Is.Null);
        }

        [Test]
        public static void TestCustomResidue()
        {
            string name = "FakeNucleotide";
            char oneLetter = 'F';
            string symbol = "Fke";
            string chemicalFormula = "C5H5N2O2";
            var fakeNucleotide = new Nucleotide(name, oneLetter, symbol, ChemicalFormula.ParseFormula(chemicalFormula));

            Nucleotide.AddResidue(name, oneLetter, symbol, chemicalFormula);

            // test new nucleotide is within dictionary
            if (Nucleotide.TryGetResidue('F', out Nucleotide outTide))
            {
                Assert.That(fakeNucleotide.Equals(outTide));
            }
            else
                Assert.Fail();

            if (Nucleotide.TryGetResidue("Fke", out outTide))
            {
                Assert.That(fakeNucleotide.Equals(outTide));
            }
            else
                Assert.Fail();

            // test false result in TryGetResidue
            if (Nucleotide.TryGetResidue('P', out outTide))
                Assert.Fail();

            if (Nucleotide.TryGetResidue("Taco", out outTide))
                Assert.Fail();
        }

        [Test]
        public void TestEquality()
        {
            Assert.That(Nucleotide.TryGetResidue('A', out Nucleotide a));
            Assert.That(Nucleotide.TryGetResidue("Ade", out Nucleotide a2));
            Assert.That(Nucleotide.TryGetResidue("U", out Nucleotide u));
            Assert.That(Nucleotide.TryGetResidue("Ura", out Nucleotide u2));

            Assert.That(a.Equals(a));
            Assert.That(a.Equals(a2));
            Assert.That(a.GetHashCode(), Is.EqualTo(a2.GetHashCode()));
            Assert.That(u.Equals(u2));
            Assert.That(u.GetHashCode(), Is.EqualTo(u2.GetHashCode()));
            Assert.That(!a.Equals(u2));
            Assert.That(a.GetHashCode(), Is.Not.EqualTo(u.GetHashCode()));
            Assert.That(!u.Equals(a2));
            Assert.That(u.GetHashCode(), Is.Not.EqualTo(a.GetHashCode()));
            Assert.That(!u.Equals(null));
            Assert.That(a.Equals((object)a2));
            Assert.That(a.Equals((object)a));
            Assert.That(u.Equals((object)u2));
            Assert.That(!a.Equals((object)u2));
            Assert.That(!u.Equals((object)a2));
            Assert.That(!u.Equals((object)null));
            Assert.That(!u.Equals((object)new Action(() => { })));
        }
    }
}
