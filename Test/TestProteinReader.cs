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
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class TestProteinReader
    {
        #region Public Methods

        [Test]
        public static void Compare_protein_properties()
        {
            DatabaseReference d = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference dd = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference de = new DatabaseReference("asdf", "asdefg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference df = new DatabaseReference("asddf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference dg = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("babb", "ccc") });
            DatabaseReference dh = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "cccf") });
            Assert.True(dd.Equals(d));
            Assert.False(de.Equals(d));
            Assert.False(df.Equals(d));
            Assert.False(dg.Equals(d));
            Assert.False(dh.Equals(d));
            Assert.AreEqual(5, new HashSet<DatabaseReference> { d, dd, de, df, dg, dh }.Count);

            SequenceVariation s = new SequenceVariation(1, "hello", "hey", "hi");
            SequenceVariation sv = new SequenceVariation(1, "hello", "hey", "hi");
            SequenceVariation sss = new SequenceVariation(2, "hallo", "hey", "hi");
            SequenceVariation ssss = new SequenceVariation(1, "hello", "heyy", "hi");
            SequenceVariation sssss = new SequenceVariation(1, "hello", "hey", "hii");
            Assert.True(s.Equals(sv));
            Assert.False(s.Equals(sss));
            Assert.False(s.Equals(ssss));
            Assert.False(s.Equals(sssss));
            Assert.AreEqual(4, new HashSet<SequenceVariation> { s, sv, sss, ssss, sssss }.Count);

            DisulfideBond b = new DisulfideBond(1, "hello");
            DisulfideBond bb = new DisulfideBond(1, "hello");
            DisulfideBond bbb = new DisulfideBond(1, 2, "hello");
            DisulfideBond bbbb = new DisulfideBond(1, 2, "hello");
            DisulfideBond ba = new DisulfideBond(1, 3, "hello");
            DisulfideBond baa = new DisulfideBond(2, 2, "hello");
            DisulfideBond baaa = new DisulfideBond(1, 2, "hallo");
            Assert.AreEqual(b, bb);
            Assert.AreEqual(bbb, bbbb);
            Assert.AreNotEqual(b, bbb);
            Assert.AreNotEqual(ba, bbb);
            Assert.AreNotEqual(baa, bbb);
            Assert.AreNotEqual(baaa, bbb);
            Assert.AreEqual(5, new HashSet<DisulfideBond> { b, bb, bbb, bbbb, ba, baa, baaa }.Count);

            ProteolysisProduct pp = new ProteolysisProduct(1, 1, "hello");
            ProteolysisProduct paaa = new ProteolysisProduct(1, 1, "hello");
            ProteolysisProduct p = new ProteolysisProduct(null, null, "hello");
            ProteolysisProduct ppp = new ProteolysisProduct(1, 2, "hello");
            ProteolysisProduct pa = new ProteolysisProduct(2, 1, "hello");
            ProteolysisProduct paa = new ProteolysisProduct(1, 1, "hallo");
            Assert.AreEqual(pp, paaa);
            Assert.AreNotEqual(p, pp);
            Assert.AreNotEqual(pp, ppp);
            Assert.AreNotEqual(pp, pa);
            Assert.AreNotEqual(pp, paa);
            Assert.AreEqual(5, new HashSet<ProteolysisProduct> { p, pp, ppp, pa, paa, paaa }.Count);
        }

        [Test]
        public static void Merge_a_couple_proteins()
        {
            ModificationMotif.TryGetMotif("A", out ModificationMotif motif);
            Protein p = new Protein(
                "ASEQUENCE",
                "id",
                isContaminant: false,
                isDecoy: false,
                name: "name",
                full_name: "full_name",
                gene_names: new List<Tuple<string, string>> { new Tuple<string, string>("gene", "name") },
                databaseReferences: new List<DatabaseReference> { new DatabaseReference("ref", "id", new List<Tuple<string, string>> { new Tuple<string, string>("type", "property") }) },
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, 2, "A", "B", "var") },
                proteolysisProducts: new List<ProteolysisProduct> { new ProteolysisProduct(1, 2, "prod") },
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { new ModificationWithMass("mod", "type", motif, TerminusLocalization.Any, 1, null, null, null) } } }
                );

            Protein p2 = new Protein(
                "ASEQUENCE",
                "id",
                isContaminant: false,
                isDecoy: false,
                name: "name",
                full_name: "full_name",
                gene_names: new List<Tuple<string, string>> { new Tuple<string, string>("gene", "name") },
                databaseReferences: new List<DatabaseReference> { new DatabaseReference("ref", "id", new List<Tuple<string, string>> { new Tuple<string, string>("type", "property") }) },
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, 2, "A", "B", "var") },
                proteolysisProducts: new List<ProteolysisProduct> { new ProteolysisProduct(1, 2, "prod") },
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { new ModificationWithMass("mod2", "type", motif, TerminusLocalization.Any, 10, null, null, null) } } }
                );

            List<Protein> merged = ProteinDbLoader.Merge_proteins(new List<Protein> { p, p2 }).ToList();
            Assert.AreEqual(1, merged.Count);
            Assert.AreEqual(1, merged.First().DatabaseReferences.Count());
            Assert.AreEqual(1, merged.First().GeneNames.Count());
            Assert.AreEqual(1, merged.First().SequenceVariations.Count());
            Assert.AreEqual(1, merged.First().ProteolysisProducts.Count());
            Assert.AreEqual(p.OneBasedPossibleLocalizedModifications.First().Value.First().GetHashCode(), p2.OneBasedPossibleLocalizedModifications.First().Value.First().GetHashCode());
            Assert.AreNotEqual(p.OneBasedPossibleLocalizedModifications.First().Value.First(), p2.OneBasedPossibleLocalizedModifications.First().Value.First());
            Assert.AreEqual(1, merged.First().OneBasedPossibleLocalizedModifications.Count());
            Assert.AreEqual(2, merged.First().OneBasedPossibleLocalizedModifications.First().Value.Count);
        }

        [Test]
        public static void XmlTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",  null, null,TerminusLocalization.Any,null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml"), true, true, nice, false, null, out Dictionary<string, Modification> un);

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
        public static void SeqVarXmlTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",  null, null,TerminusLocalization.Any,null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"seqvartests.xml"), true, true, nice, false, null, out Dictionary<string, Modification> un);

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
        public static void DisulfideXmlTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",  null, null,TerminusLocalization.Any,null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"disulfidetests.xml"), true, true, nice, false, null, out Dictionary<string, Modification> un);

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
        public static void XmlTest_2entry()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",  null, null,TerminusLocalization.Any,null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), true, true, nice, false, null, out Dictionary<string, Modification> un);

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
        public static void XmlGzTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",  null, null,TerminusLocalization.Any,null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml.gz"), true, true, nice, false, null, out Dictionary<string, Modification> un);

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
        public static void XmlFunkySequenceTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",  null, null,TerminusLocalization.Any,null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"fake_h4.xml"), true, true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('S', ok[0][0]);
            Assert.AreEqual('G', ok[1][0]);
        }

        [Test]
        public static void XmlModifiedStartTest()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",  null, null,TerminusLocalization.Any,null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"modified_start.xml"), true, true, nice, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('M', ok[0][0]);
            Assert.AreEqual('M', ok[1][0]);
            Assert.AreEqual(1, ok[1].OneBasedPossibleLocalizedModifications[1].Count);
        }

        [Test]
        public static void FastaTest()
        {
            List<Protein> prots = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"fasta.fasta"), true, true, false, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_gene_expression);
            Assert.AreEqual("P62805", prots.First().Accession);
            Assert.AreEqual("H4_HUMAN Histone H4", prots.First().FullName);
            Assert.AreEqual("HIST1H4A", prots.First().GeneNames.First().Item2);
        }

        [Test]
        public static void Load_fasta_handle_tooHigh_indices()
        {
            ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"bad.fasta"), true, true, false, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_fullName_expression, ProteinDbLoader.uniprot_accession_expression, ProteinDbLoader.uniprot_gene_expression);
        }

        [Test]
        public static void Read_xml_mod_collision()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("N-acetylserine", "one",  null, TerminusLocalization.Any, null),
                new ModificationWithLocation("N-acetylserine", "two",  null, TerminusLocalization.Any, null)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml"), true, true, nice, false, new List<string>(), out Dictionary<string, Modification> un);
            Assert.True(ok[0].OneBasedPossibleLocalizedModifications.Any(kv => kv.Value.Count > 1));
            Assert.True(ok[0].OneBasedPossibleLocalizedModifications[2].Select(m => m.id).Contains("N-acetylserine"));
        }

        [Test]
        public static void Read_xml_exclude_mods()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("N-acetylserine", "exclude_me",  null, TerminusLocalization.Any, null)
            };

            var ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml.xml"), true, true, nice, false, new string[] { "exclude_me" }, out Dictionary<string, Modification> un);
            Assert.False(ok2[0].OneBasedPossibleLocalizedModifications[2].Select(m => m.id).Contains("N-acetylserine"));
        }

        [Test]
        public static void CompareOxidationWithAndWithoutCf()
        {
            string aString =
@"ID   Methionine (R)-sulfoxide
AC   PTM-0480
FT   MOD_RES
TG   Methionine.
PA   Amino acid side chain.
PP   Anywhere.
CF   O1
MM   15.994915
MA   16.00
LC   Intracellular localisation.
TR   Eukaryota; taxId:2759 (Eukaryota).
KW   Oxidation.
DR   RESID; AA0581.
DR   PSI-MOD; MOD:00720.
//";
            var a = PtmListLoader.ReadModsFromString(aString).First();

            string bString =
@"ID   Oxidation of M
TG   M
PP   Anywhere.
MT   Common Variable
CF   O1
//";
            var b = PtmListLoader.ReadModsFromString(bString).First();

            Assert.IsTrue(Math.Abs((a as ModificationWithMass).monoisotopicMass - (b as ModificationWithMass).monoisotopicMass) < 1e-6);
            Assert.IsTrue(Math.Abs((a as ModificationWithMass).monoisotopicMass - (b as ModificationWithMass).monoisotopicMass) > 1e-7);
        }

        [Test]
        public static void TestKeywordAugmentation()
        {
            string bString =
@"ID   Oxidation
TG   M or R
PP   Anywhere.
MT   Common Variable
CF   O1
//";
            var a = PtmListLoader.ReadModsFromString(bString).First();
            var b = PtmListLoader.ReadModsFromString(bString).Last();

            Assert.AreEqual("Oxidation on M", a.id);
            Assert.AreEqual("Oxidation", (a as ModificationWithMass).keywords.First());
            Assert.AreEqual("Oxidation on R", b.id);
            Assert.AreEqual("Oxidation", (b as ModificationWithMass).keywords.First());
        }

        #endregion Public Methods
    }
}