// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (TestRange.cs) is part of MassSpectrometry.Tests.
//
// MassSpectrometry.Tests is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry.Tests is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry.Tests. If not, see <http://www.gnu.org/licenses/>.

using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public sealed class TestProteinReader
    {

        #region Public Methods

        [Test]
        public void XmlTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml"), true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('M', ok[0][0]);
            Assert.AreEqual('M', ok[1][0]);

            Assert.AreEqual("P62805|H4_HUMAN|Histone H4", ok[0].FullDescription);
            Assert.AreEqual("DECOY_P62805|H4_HUMAN|Histone H4", ok[1].FullDescription);
            Assert.AreEqual("ENST00000244537", ok[0].DatabaseReferences.First(dbRef => dbRef.Type == "Ensembl").Id);
            Assert.AreEqual("protein sequence ID", ok[0].DatabaseReferences.First(dbRef => dbRef.Type == "Ensembl").Properties.First().Item1);
            Assert.AreEqual("ENSP00000244537", ok[0].DatabaseReferences.First(dbRef => dbRef.Type == "Ensembl").Properties.First().Item2);
            Assert.AreEqual(42, ok[0].GeneNames.Count());
            Assert.AreEqual(14, ok[0].GeneNames.Where(t => t.Item1 == "primary").Count());
            Assert.AreEqual("HIST1H4A", ok[0].GeneNames.Where(t => t.Item1 == "primary").First().Item2);
            Assert.AreEqual(23, ok[0].DatabaseReferences.Count(dbRef => dbRef.Type == "Ensembl"));
            Assert.AreEqual(0, ok[0].DisulfideBonds.Count());
            Assert.AreEqual(1, ok[0].SequenceVariations.Count());
            Assert.AreEqual(1, ok[1].SequenceVariations.Count()); // decoys get the same sequence variations
            Assert.AreEqual(64, ok[0].SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(64, ok[0].SequenceVariations.First().OneBasedEndPosition);
            Assert.AreEqual(103 - 64 + 2, ok[1].SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(103 - 64 + 2, ok[1].SequenceVariations.First().OneBasedEndPosition);
            Assert.AreNotEqual(ok[0].SequenceVariations.First().Description, ok[1].SequenceVariations.First().Description); //decoys and target variations don't have the same desc.
        }

        [Test]
        public void SeqVarXmlTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"seqvartests.xml"), true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('M', ok[0][0]);
            Assert.AreEqual('M', ok[1][0]);
            List<SequenceVariation> seqvar0 = ok[0].SequenceVariations.ToList();
            List<SequenceVariation> seqvar1 = ok[1].SequenceVariations.ToList();
            Assert.AreEqual(seqvar0.Count + 1, seqvar1.Count);
            Assert.AreEqual('M', ok[0].SequenceVariations.First().OriginalSequence[0]);
            Assert.AreEqual('M', ok[0].SequenceVariations.First().VariantSequence[0]);
            Assert.AreEqual('A', ok[1].SequenceVariations.First().OriginalSequence[0]);
            Assert.AreEqual('P', ok[1].SequenceVariations.First().VariantSequence[0]);
            Assert.AreEqual('M', seqvar0[1].OriginalSequence[0]);
            Assert.AreEqual("", seqvar1[1].VariantSequence);
            foreach (SequenceVariation s in seqvar0)
            {
                Assert.AreEqual(s.OriginalSequence, ok[0].BaseSequence.Substring(s.OneBasedBeginPosition - 1, s.OneBasedEndPosition - s.OneBasedBeginPosition + 1));
            }
            foreach (SequenceVariation s in seqvar1)
            {
                Assert.AreEqual(s.OriginalSequence, ok[1].BaseSequence.Substring(s.OneBasedBeginPosition - 1, s.OneBasedEndPosition - s.OneBasedBeginPosition + 1));
            }
            Assert.AreNotEqual(ok[0].SequenceVariations.First().Description, ok[1].SequenceVariations.First().Description); //decoys and target variations don't have the same desc.
        }

        [Test]
        public void DisulfideXmlTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"disulfidetests.xml"), true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('M', ok[0][0]);
            Assert.AreEqual('M', ok[1][0]);

            Assert.AreEqual(3, ok[0].DisulfideBonds.Count());
            Assert.AreEqual(31, ok[0].DisulfideBonds.First().OneBasedBeginPosition);
            Assert.AreEqual(94, ok[0].DisulfideBonds.First().OneBasedEndPosition);
            Assert.AreEqual(93, ok[0].DisulfideBonds.ElementAt(2).OneBasedBeginPosition);
            Assert.AreEqual(93, ok[0].DisulfideBonds.ElementAt(2).OneBasedEndPosition);

            Assert.AreEqual(3, ok[1].DisulfideBonds.Count());
            Assert.AreEqual(78, ok[1].DisulfideBonds.First().OneBasedBeginPosition);
            Assert.AreEqual(15, ok[1].DisulfideBonds.First().OneBasedEndPosition);
            Assert.AreEqual(16, ok[1].DisulfideBonds.ElementAt(2).OneBasedBeginPosition);
            Assert.AreEqual(16, ok[1].DisulfideBonds.ElementAt(2).OneBasedEndPosition);
            Assert.AreNotEqual(ok[0].DisulfideBonds.First().Description, ok[1].DisulfideBonds.First().Description); //decoys and target disulfide bonds don't have the same desc.
        }

        [Test]
        public void XmlTest_2entry()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.True(ok.All(p => p.ProteolysisProducts.All(d => d.OneBasedBeginPosition == null || d.OneBasedBeginPosition > 0)));

            Assert.True(ok.All(p => p.ProteolysisProducts.All(d => d.OneBasedEndPosition == null || d.OneBasedEndPosition <= p.Length)));

            Assert.False(ok.All(p => p.BaseSequence.Contains(" ")));
            Assert.False(ok.All(p => p.BaseSequence.Contains("\t")));
            Assert.False(ok.All(p => p.BaseSequence.Contains("\n")));

            //GoTerm checks
            List<Protein> targets = ok.Where(p => !p.IsDecoy).ToList();
            Assert.AreEqual(2, targets.Count);
            Assert.AreEqual(1, targets[0].DatabaseReferences.Count(dbRef => dbRef.Type == "EnsemblFungi"));
            Assert.AreEqual(1, targets[1].DatabaseReferences.Count(dbRef => dbRef.Type == "EnsemblFungi"));
        }

        [Test]
        public void XmlGzTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml.gz"), true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('M', ok[0][0]);
            Assert.AreEqual('M', ok[1][0]);

            Assert.AreEqual("P62805|H4_HUMAN|Histone H4", ok[0].FullDescription);
            Assert.AreEqual("DECOY_P62805|H4_HUMAN|Histone H4", ok[1].FullDescription);
            Assert.AreEqual("ENST00000244537", ok[0].DatabaseReferences.First(dbRef => dbRef.Type == "Ensembl").Id);
            Assert.AreEqual("protein sequence ID", ok[0].DatabaseReferences.First(dbRef => dbRef.Type == "Ensembl").Properties.First().Item1);
            Assert.AreEqual("ENSP00000244537", ok[0].DatabaseReferences.First(dbRef => dbRef.Type == "Ensembl").Properties.First().Item2);
            Assert.AreEqual(42, ok[0].GeneNames.Count());
            Assert.AreEqual(14, ok[0].GeneNames.Where(t => t.Item1 == "primary").Count());
            Assert.AreEqual("HIST1H4A", ok[0].GeneNames.Where(t => t.Item1 == "primary").First().Item2);
            Assert.AreEqual(23, ok[0].DatabaseReferences.Count(dbRef => dbRef.Type == "Ensembl"));
        }

        [Test]
        public void XmlFunkySequenceTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"fake_h4.xml"), true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('S', ok[0][0]);
            Assert.AreEqual('G', ok[1][0]);
        }

        [Test]
        public void XmlModifiedStartTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"modified_start.xml"), true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('M', ok[0][0]);
            Assert.AreEqual('M', ok[1][0]);
            Assert.AreEqual(1, ok[1].OneBasedPossibleLocalizedModifications[1].Count);
        }

        [Test]
        public void FastaTest()
        {
            List<Protein> prots = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"fasta.fasta"), true, false, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_gene_expression);
            Assert.AreEqual("P62805", prots.First().Accession);
            Assert.AreEqual("H4_HUMAN Histone H4", prots.First().FullName);
            Assert.AreEqual("HIST1H4A", prots.First().GeneNames.First().Item2);
        }

        [Test]
        public void Load_fasta_handle_tooHigh_indices()
        {
            ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"bad.fasta"), true, false, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_gene_expression);
        }

        [Test]
        public void Read_xml_mod_collision()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("N-acetylserine", null, null, ModificationSites.S, null, "one"),
                new ModificationWithLocation("N-acetylserine", null, null, ModificationSites.S, null, "two")
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml"), true, nice, false, new List<string>(), out Dictionary<string, Modification> un);
            Assert.True(ok[0].OneBasedPossibleLocalizedModifications.Any(kv => kv.Value.Count > 1));
            Assert.True(ok[0].OneBasedPossibleLocalizedModifications[2].Select(m => m.id).Contains("N-acetylserine"));
        }

        [Test]
        public void Read_xml_exclude_mods()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("N-acetylserine", null, null, ModificationSites.S, null, "exclude_me")
            };

            var ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml"), true, nice, false, new string[] { "exclude_me" }, out Dictionary<string, Modification> un);
            Assert.False(ok2[0].OneBasedPossibleLocalizedModifications[2].Select(m => m.id).Contains("N-acetylserine"));
        }

        #endregion Public Methods

    }
}