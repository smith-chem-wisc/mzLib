using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestCompactPeptide
    {
        [Test]
        public static void TestCompactPeptideMasses_UnmodifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);

            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);
            Assert.AreEqual(2, aCompactPeptide.NTerminalMasses.Length);
            Assert.AreEqual(2, aCompactPeptide.CTerminalMasses.Length);

            //evaluate N-terminal masses
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 97, 226 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 101, 230 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);

        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);

            //evaluate N-terminal masses
            Assert.AreEqual(2, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 177, 306 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(2, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 101, 230 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

        [Test]
        public static void TestCompactPeptideMasses_cTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);

            //evaluate N-terminal masses
            Assert.AreEqual(2, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 97, 226 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(2, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 181, 310 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

        [Test]
        public static void TestCompactPeptideMasses_internallyModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);

            //evaluate N-terminal masses
            Assert.AreEqual(2, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 97, 306 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(2, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 101, 310 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);

            //evaluate N-terminal masses
            Assert.AreEqual(4, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 177, 306, 79, 208 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(2, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 101, 230 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

        [Test]
        public static void TestCompactPeptideMasses_cTerminalModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);

            //evaluate N-terminal masses
            Assert.AreEqual(2, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 97, 226 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(4, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 181, 310, 83, 212 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

        [Test]
        public static void TestCompactPeptideMasses_internallyModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);

            //evaluate N-terminal masses
            Assert.AreEqual(3, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 97, 306, 208 }; 
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(4, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 101, 310, 212 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide_NeutralLoss_DissociationTypes_AnyActivationType_and_HCD()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);//this dissociation type was the one used in the experiment

            //evaluate N-terminal masses
            Assert.AreEqual(4, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 177, 306, 79, 208 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(2, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 101, 230 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide_NeutralLoss_DissociationTypes_CID_and_HCD()//there should be no added neutral losses in this case becuase the allowed dissociation type doesn't match the dissociation type used in the experiment
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_id: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.HCD);//this dissociation type was the one used in the experiment

            //evaluate N-terminal masses
            Assert.AreEqual(2, aCompactPeptide.NTerminalMasses.Length);
            List<int> nTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.NTerminalMasses)
            {
                nTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_n_Values = new List<int> { 177, 306 };
            Assert.IsTrue(nTerminalMassListRounded.Except(actual_n_Values).ToList().Count == 0);

            //evaluate C-terminal masses
            Assert.AreEqual(2, aCompactPeptide.CTerminalMasses.Length);
            List<int> cTerminalMassListRounded = new List<int>();
            foreach (double mass in aCompactPeptide.CTerminalMasses)
            {
                cTerminalMassListRounded.Add((int)Math.Round(mass, MidpointRounding.AwayFromZero));
            }
            List<int> actual_c_Values = new List<int> { 101, 230 };
            Assert.IsTrue(cTerminalMassListRounded.Except(actual_c_Values).ToList().Count == 0);
        }

    }
}