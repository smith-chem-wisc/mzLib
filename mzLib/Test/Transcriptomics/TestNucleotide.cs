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
    internal class TestNucleotide
    {
        internal record NucleotideTestCase(Nucleotide Nucleotide, string Name, char OneLetterCode, string Symbol, ChemicalFormula Formula, double Mass,
            ChemicalFormula nucleosideFormula);

        internal static IEnumerable<NucleotideTestCase> GetNucleotideTestCases()
        {
            Loaders.LoadElements();

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
