using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Digestion hands every modified form of a digestion product the product's own base sequence string and residue
    /// mass sum, instead of each form taking its own substring and summing its own residues. These tests hold each form
    /// to what the previous code computed for it, kept verbatim below: the same base sequence, and the same monoisotopic
    /// mass bit for bit - and to a form built through the public constructor, which still computes everything itself.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class DigestionProductInvariantsTests
    {
        [Test]
        public static void EveryFormHasTheBaseSequenceAndMassItComputedBefore()
        {
            var random = new Random(20261015);
            long forms = 0, products = 0, sharedForms = 0;

            for (int trial = 0; trial < 400; trial++)
            {
                DigestionParams digestionParams = RandomDigestionInputs.RandomDigestionParams(random);
                (List<Modification> fixedMods, List<Modification> variableMods) = RandomDigestionInputs.RandomSearchMods(random);
                Protein protein = RandomDigestionInputs.RandomProtein(random, trial, 5, 400, isDecoy: random.Next(4) == 0);

                List<PeptideWithSetModifications> digested = protein.Digest(digestionParams, fixedMods, variableMods).ToList();
                PeptideWithSetModifications previous = null;
                foreach (PeptideWithSetModifications form in digested)
                {
                    string expectedBase = BaseSequenceBefore(form);
                    Assert.That(form.BaseSequence, Is.EqualTo(expectedBase));
                    Assert.That(BitConverter.DoubleToInt64Bits(form.MonoisotopicMass),
                        Is.EqualTo(BitConverter.DoubleToInt64Bits(MonoisotopicMassBefore(expectedBase, form.AllModsOneIsNterminus))),
                        $"mass of {form.FullSequence}");

                    // the same form built the public way, which takes its own substring and sums its own residues
                    var reference = new PeptideWithSetModifications(protein, digestionParams, form.OneBasedStartResidue, form.OneBasedEndResidue,
                        form.CleavageSpecificityForFdrCategory, form.PeptideDescription, form.MissedCleavages,
                        new Dictionary<int, Modification>(form.AllModsOneIsNterminus), form.NumFixedMods);
                    Assert.That(form.BaseSequence, Is.EqualTo(reference.BaseSequence));
                    Assert.That(form.FullSequence, Is.EqualTo(reference.FullSequence));
                    Assert.That(BitConverter.DoubleToInt64Bits(form.MonoisotopicMass), Is.EqualTo(BitConverter.DoubleToInt64Bits(reference.MonoisotopicMass)));
                    Assert.That(form, Is.EqualTo(reference));

                    // forms of one product follow one another, and now hold the very same base sequence string
                    bool sameProduct = previous != null && previous.OneBasedStartResidue == form.OneBasedStartResidue
                        && previous.OneBasedEndResidue == form.OneBasedEndResidue && previous.CleavageSpecificityForFdrCategory == form.CleavageSpecificityForFdrCategory;
                    if (sameProduct && ReferenceEquals(previous.BaseSequence, form.BaseSequence))
                    {
                        sharedForms++;
                    }
                    else if (!sameProduct)
                    {
                        products++;
                    }
                    previous = form;
                    forms++;
                }
            }

            Assert.That(forms, Is.GreaterThan(20000));
            Assert.That(forms - products, Is.GreaterThan(5000), "products must have several forms for sharing to be exercised");
            Assert.That(sharedForms, Is.GreaterThan(5000));
        }

        [Test]
        public static void UnknownResiduesStillGiveAnUnknownMass()
        {
            var protein = new Protein("PEPXIDEK", "P1");
            var digestionParams = new DigestionParams("trypsin", minPeptideLength: 1);
            PeptideWithSetModifications form = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()).Single();
            Assert.That(double.IsNaN(form.MonoisotopicMass));
            Assert.That(BitConverter.DoubleToInt64Bits(form.MonoisotopicMass),
                Is.EqualTo(BitConverter.DoubleToInt64Bits(MonoisotopicMassBefore(form.BaseSequence, form.AllModsOneIsNterminus))));
        }

        #region the previous code, verbatim

        // DigestionProduct.BaseSequence, for a form constructed without a base sequence
        private static string BaseSequenceBefore(PeptideWithSetModifications form) =>
            form.Parent.BaseSequence.Substring(form.OneBasedStartResidue - 1, form.OneBasedEndResidue - form.OneBasedStartResidue + 1);

        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        // PeptideWithSetModifications.MonoisotopicMass
        private static double MonoisotopicMassBefore(string BaseSequence, Dictionary<int, Modification> AllModsOneIsNterminus)
        {
            double monoMass = WaterMonoisotopicMass;

            foreach (var mod in AllModsOneIsNterminus.Values)
            {
                monoMass += mod.MonoisotopicMass.Value;
            }
            monoMass += BaseSequence.Sum(b => Residue.ResidueMonoisotopicMass[b]);

            return (double)ClassExtensions.RoundedDouble(monoMass);
        }

        #endregion
    }
}
