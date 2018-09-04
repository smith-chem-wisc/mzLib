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

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);
            
            //evaluate N-terminal masses
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 1))));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 101, 230, 327 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 1))));
        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);
            
            //evaluate N-terminal masses
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 101, 230, 407 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
        }

        [Test]
        public static void TestCompactPeptideMasses_cTerminalModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);
            
            //evaluate N-terminal masses
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 1))));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 181, 310, 407 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 1))));
        }

        [Test]
        public static void TestCompactPeptideMasses_internallyModifiedPeptide()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"));
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 306 };
            HashSet<int> foundNTerminalMasses = new HashSet<int>(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());

            Assert.That(expectedNTerminalMasses.SetEquals(foundNTerminalMasses));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 101, 310, 407 };
            HashSet<int> foundCTerminalMasses = new HashSet<int>(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());

            Assert.That(expectedCTerminalMasses.SetEquals(foundCTerminalMasses));
        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 101, 230, 407 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
        }

        [Test]
        public static void TestCompactPeptideMasses_cTerminalModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);

            //var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> foundNTerminalMasses = new HashSet<int>(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 226 };

            Assert.That(expectedNTerminalMasses.SetEquals(foundNTerminalMasses));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> foundCTerminalMasses = new HashSet<int>(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0)).ToList());
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 181, 310, 407 };

            Assert.That(expectedCTerminalMasses.SetEquals(foundCTerminalMasses));
        }

        [Test]
        public static void TestCompactPeptideMasses_internallyModifiedPeptide_NeutralLoss()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.HCD, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);

            var allFragmentNeutralMasses = aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both);
            
            //evaluate N-terminal masses
            var n = allFragmentNeutralMasses.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.N).ToList();
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 97, 306, 208 };
            Assert.That(expectedNTerminalMasses.SetEquals(n.Select(v => (int)Math.Round(v.NeutralMass, 1))));

            //evaluate C-terminal masses
            var c = allFragmentNeutralMasses.Where(f => f.TerminusFragment.Terminus == FragmentationTerminus.C).ToList();
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 119, 328, 230 };
            Assert.That(expectedCTerminalMasses.SetEquals(c.Select(v => (int)Math.Round(v.NeutralMass, 1))));
        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide_NeutralLoss_DissociationTypes_AnyActivationType_and_HCD()
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 101, 230, 407 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 0))));
        }

        [Test]
        public static void TestCompactPeptideMasses_nTerminalModifiedPeptide_NeutralLoss_DissociationTypes_CID_and_HCD()//there should be no added neutral losses in this case becuase the allowed dissociation type doesn't match the dissociation type used in the experiment
        {
            Protein p = new Protein("PET", "accession");
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification phosphorylation = new Modification(_originalId: "phospho", _modificationType: "CommonBiological", _target: motif, _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.CID, new List<double> { 0, ChemicalFormula.ParseFormula("H3O4P1").MonoisotopicMass } } });
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);
            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification> { phosphorylation }, new List<Modification>()).First();

            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both);

            //evaluate N-terminal masses
            var nTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.N);
            HashSet<int> expectedNTerminalMasses = new HashSet<int> { 177, 306 };
            Assert.That(expectedNTerminalMasses.SetEquals(nTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 1))));

            //evaluate C-terminal masses
            var cTerminalMasses = aCompactPeptide.TerminalMasses.Where(v => v.Terminus == FragmentationTerminus.C);
            HashSet<int> expectedCTerminalMasses = new HashSet<int> { 101, 230, 407 };
            Assert.That(expectedCTerminalMasses.SetEquals(cTerminalMasses.Select(v => (int)Math.Round(v.NeutralMass, 1))));
        }
    }
}