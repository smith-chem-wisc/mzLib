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
        private static List<Modification> UniProtPtms;

        [OneTimeSetUp]
        public static void SetUpModifications()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();
        }

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
                fullName: "full_name",
                geneNames: new List<Tuple<string, string>> { new Tuple<string, string>("gene", "name") },
                databaseReferences: new List<DatabaseReference> { new DatabaseReference("ref", "id", new List<Tuple<string, string>> { new Tuple<string, string>("type", "property") }) },
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, 2, "A", "B", "var") },
                proteolysisProducts: new List<ProteolysisProduct> { new ProteolysisProduct(1, 2, "prod") },
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { new Modification("mod", null, "type", null, motif, "Anywhere.", null, 1, null, null, null, null, null, null) } } }
                );

            Protein p2 = new Protein(
                "ASEQUENCE",
                "id",
                isContaminant: false,
                isDecoy: false,
                name: "name",
                fullName: "full_name",
                geneNames: new List<Tuple<string, string>> { new Tuple<string, string>("gene", "name") },
                databaseReferences: new List<DatabaseReference> { new DatabaseReference("ref", "id", new List<Tuple<string, string>> { new Tuple<string, string>("type", "property") }) },
                sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, 2, "A", "B", "var") },
                proteolysisProducts: new List<ProteolysisProduct> { new ProteolysisProduct(1, 2, "prod") },
                oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { new Modification("mod", null, "type", null, motif, "Anywhere.", null, 10, null, null, null, null, null, null) } } }
                );

            List<Protein> merged = ProteinDbLoader.MergeProteins(new List<Protein> { p, p2 }).ToList();
            Assert.AreEqual(1, merged.Count);
            Assert.AreEqual(1, merged.First().DatabaseReferences.Count());
            Assert.AreEqual(1, merged.First().GeneNames.Count());
            Assert.AreEqual(1, merged.First().SequenceVariations.Count());
            Assert.AreEqual(1, merged.First().ProteolysisProducts.Count());
            Assert.AreNotEqual(p.OneBasedPossibleLocalizedModifications.First().Value.First(), p2.OneBasedPossibleLocalizedModifications.First().Value.First());
            Assert.AreEqual(1, merged.First().OneBasedPossibleLocalizedModifications.Count());
            Assert.AreEqual(2, merged.First().OneBasedPossibleLocalizedModifications.First().Value.Count);
        }

        [Test]
        public static void XmlTest()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml"), 
                true, DecoyType.Reverse, UniProtPtms, false, null, out var un);

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
            Assert.AreEqual("Homo sapiens", ok[1].Organism);
        }

        [Test]
        public static void SeqVarXmlTest()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "seqvartests.xml"),
                true, DecoyType.Reverse, UniProtPtms, false, null, out var un);

            var target = ok.First(p => !p.IsDecoy);
            Protein decoy = ok.Where(p => p.IsDecoy && p.SequenceVariations.Count() > 0).First();
            
            Assert.AreEqual('M', target[0]);
            Assert.AreEqual('M', decoy[0]);
            List<SequenceVariation> seqvar0 = target.SequenceVariations.ToList();
            List<SequenceVariation> seqvar1 = decoy.SequenceVariations.ToList();
            Assert.AreEqual(seqvar0.Count + 1, seqvar1.Count);
            Assert.AreEqual('M', target.SequenceVariations.First().OriginalSequence[0]);
            Assert.AreEqual('M', target.SequenceVariations.First().VariantSequence[0]);
            Assert.AreEqual('A', decoy.SequenceVariations.First().OriginalSequence[0]);
            Assert.AreEqual('P', decoy.SequenceVariations.First().VariantSequence[0]);
            Assert.AreEqual('M', seqvar0[1].OriginalSequence[0]);
            Assert.AreEqual("", seqvar1[1].VariantSequence);
            foreach (SequenceVariation s in seqvar0)
            {
                Assert.AreEqual(s.OriginalSequence, target.BaseSequence.Substring(s.OneBasedBeginPosition - 1, s.OneBasedEndPosition - s.OneBasedBeginPosition + 1));
            }
            foreach (SequenceVariation s in seqvar1)
            {
                Assert.AreEqual(s.OriginalSequence, decoy.BaseSequence.Substring(s.OneBasedBeginPosition - 1, s.OneBasedEndPosition - s.OneBasedBeginPosition + 1));
            }
            Assert.AreNotEqual(target.SequenceVariations.First().Description, decoy.SequenceVariations.First().Description); //decoys and target variations don't have the same desc.
        }

        [Test]
        public static void DisulfideXmlTest()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"disulfidetests.xml"), 
                true, DecoyType.Reverse, UniProtPtms, false, null, out Dictionary<string, Modification> un);

            Assert.AreEqual('M', ok[0][0]);
            Assert.AreEqual('M', ok[1][0]);

            Assert.AreEqual(3, ok[0].DisulfideBonds.Count());
            Assert.AreEqual(31, ok[0].DisulfideBonds.First().OneBasedBeginPosition);
            Assert.AreEqual(94, ok[0].DisulfideBonds.First().OneBasedEndPosition);
            Assert.AreEqual(93, ok[0].DisulfideBonds.ElementAt(2).OneBasedBeginPosition);
            Assert.AreEqual(93, ok[0].DisulfideBonds.ElementAt(2).OneBasedEndPosition);

            Assert.AreEqual(3, ok[1].DisulfideBonds.Count());
            Assert.AreEqual(79, ok[1].DisulfideBonds.First().OneBasedBeginPosition);
            Assert.AreEqual(16, ok[1].DisulfideBonds.First().OneBasedEndPosition);
            Assert.AreEqual(17, ok[1].DisulfideBonds.ElementAt(2).OneBasedBeginPosition);
            Assert.AreEqual(17, ok[1].DisulfideBonds.ElementAt(2).OneBasedEndPosition);
            Assert.AreNotEqual(ok[0].DisulfideBonds.First().Description, ok[1].DisulfideBonds.First().Description); //decoys and target disulfide bonds don't have the same desc.
        }

        [Test]
        public static void XmlTest_2entry()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml2.xml"),
                true, DecoyType.Reverse, UniProtPtms, false, null, out var un);

            // proteolysis products check
            Assert.True(ok.All(p => p.ProteolysisProducts.All(d => d.OneBasedBeginPosition == null || d.OneBasedBeginPosition > 0)));
            Assert.True(ok.All(p => p.ProteolysisProducts.All(d => d.OneBasedEndPosition == null || d.OneBasedEndPosition <= p.Length)));

            // base sequence check
            Assert.False(ok.All(p => p.BaseSequence.Contains(" ")));
            Assert.False(ok.All(p => p.BaseSequence.Contains("\t")));
            Assert.False(ok.All(p => p.BaseSequence.Contains("\n")));

            // GoTerm checks
            List<Protein> targets = ok.Where(p => !p.IsDecoy).ToList();
            Assert.AreEqual(2, targets.Count);
            Assert.AreEqual(1, targets[0].DatabaseReferences.Count(dbRef => dbRef.Type == "EnsemblFungi"));
            Assert.AreEqual(1, targets[1].DatabaseReferences.Count(dbRef => dbRef.Type == "EnsemblFungi"));
        }

        [Test]
        public static void XmlGzTest()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml.gz"), 
                true, DecoyType.Reverse, UniProtPtms, false, null, out var un);

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
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"fake_h4.xml"), 
                true, DecoyType.Reverse, UniProtPtms, false, null, out var un);

            Assert.AreEqual("S", ok[0].BaseSequence.Substring(0, 1));
            Assert.AreEqual("G", ok[1].BaseSequence.Substring(0, 1));

            Assert.AreEqual('S', ok[0][0]);
            Assert.AreEqual('G', ok[1][0]);
        }

        [Test]
        public static void XmlModifiedStartTest()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"modified_start.xml"), 
                true, DecoyType.Reverse, UniProtPtms, false, null, out var un);

            Assert.AreEqual("M", ok[0].BaseSequence.Substring(0, 1)); //the original protein sequence in the original order starts with 'M'
            Assert.AreEqual("M", ok[1].BaseSequence.Substring(0, 1)); //the decoy protein sequence in the reverse order from the original still starts with 'M'
            Assert.AreEqual(1, ok[1].OneBasedPossibleLocalizedModifications[1].Count); //the initial methionine of the decoy still has the mod that it's supposed to have.
        }

        [Test]
        public static void FastaTest()
        {
            List<Protein> prots = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"fasta.fasta"), true, DecoyType.Reverse, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, out var a);

            Assert.AreEqual("P62805", prots.First().Accession);
            Assert.AreEqual("H4_HUMAN", prots.First().Name);
            Assert.AreEqual("Histone H4", prots.First().FullName);
            Assert.AreEqual("HIST1H4A", prots.First().GeneNames.First().Item2);
            Assert.AreEqual("Homo sapiens", prots.First().Organism);
        }

        [Test]
        public static void BadFastaTest()
        {
            ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bad4.fasta"), true, DecoyType.Reverse, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, out var a);
            Assert.AreEqual(1, a.Count);
            ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bad3.fasta"), true, DecoyType.Reverse, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, out var b);
            Assert.AreEqual(2, b.Count);
            ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"blank.fasta"), true, DecoyType.Reverse, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, out var c);
            Assert.AreEqual(1, c.Count);
        }

        [Test]
        public static void Load_fasta_handle_tooHigh_indices()
        {
            var p = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"bad.fasta"), true, DecoyType.Reverse, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, out var a);
        }

        [Test]
        public static void Read_xml_mod_collision()
        {
            ModificationMotif.TryGetMotif("S", out ModificationMotif motif);
            var nice = new List<Modification>
            {
                new Modification(_originalId: "N-acetylserine", _modificationType: "one", _target:  motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 1),
                new Modification(_originalId: "N-acetylserine", _modificationType: "two", _target:  motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 1)
            };

            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml"), true, DecoyType.Reverse, UniProtPtms.Concat(nice), false,
                new List<string>(), out Dictionary<string, Modification> un);

            Assert.True(ok[0].OneBasedPossibleLocalizedModifications.Any(kv => kv.Value.Count > 1));

            List<string> myOriginalIds = ok[0].OneBasedPossibleLocalizedModifications[2].Select(i => i.OriginalId).ToList();

            Assert.True(myOriginalIds.Contains("N-acetylserine"));
        }

        [Test]
        [TestCase("exclude_me", false)]//the first part is the test case, the latter part is ther result of the assertion
        //[TestCase("exclude_me_not", true)]
        public static void Read_xml_exclude_mods(string excludeString, bool isExcluded)
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);


            var nice = new List<Modification>
            {
                new Modification("N-acetylserine", null, "exclude_me", null, motif, "Anywhere.", null, 10, null, null, null, null, null, null),
                new Modification("N-acetylserine", null, "exclude_me_not", null, motif, "Anywhere.", null, 10, null, null, null, null, null, null)
            };

            Assert.That(nice[0].ValidModification);

            var ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"xml.xml"), true, DecoyType.Reverse, nice, false,
                new[] { excludeString }, out Dictionary<string, Modification> un);

            List<string> modTypes = new List<string>();
            foreach (KeyValuePair<int,List<Modification>> entry in ok2[0].OneBasedPossibleLocalizedModifications)
            {
                modTypes.AddRange(entry.Value.Select(m => m.ModificationType).ToList().Distinct());
            }
            Assert.AreEqual(isExcluded, modTypes.Contains("exclude_me"));
            Assert.AreEqual(!isExcluded, modTypes.Contains("exclude_me_not"));
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
            var a = PtmListLoader.ReadModsFromString(aString, out var errorsA).First();

            string bString =
@"ID   Oxidation of M
TG   M
PP   Anywhere.
MT   Common Variable
CF   O1
//";
            var b = PtmListLoader.ReadModsFromString(bString, out var errorsB).First();

            Assert.IsTrue(Math.Abs((double)(a as Modification).MonoisotopicMass - (double)(b as Modification).MonoisotopicMass) < 1e-6);
            Assert.IsTrue(Math.Abs((double)(a as Modification).MonoisotopicMass - (double)(b as Modification).MonoisotopicMass) > 1e-7);
        }

        [Test]
        public static void TestReverseDecoyXML()
        {
            var nice = new List<Modification>();
            var ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"disulfidetests.xml"), true, DecoyType.Reverse, nice, false,
                new string[] { "exclude_me" }, out Dictionary<string, Modification> un);

            Assert.AreEqual("MALLVHFLPLLALLALWEPKPTQAFVKQHLCGPHLVEALYLVCGERGFFYTPKSRREVEDPQVEQLELGGSPGDLQTLALEVARQKRGIVDQCCTSICSLYQLENYCN", ok2[0].BaseSequence);
            Assert.AreEqual("MNCYNELQYLSCISTCCQDVIGRKQRAVELALTQLDGPSGGLELQEVQPDEVERRSKPTYFFGREGCVLYLAEVLHPGCLHQKVFAQTPKPEWLALLALLPLFHVLLA", ok2[1].BaseSequence);
            Assert.AreEqual(ok2[0].DisulfideBonds.Count(), ok2[1].DisulfideBonds.Count());
            Assert.AreEqual(ok2[0].ProteolysisProducts.Count(), ok2[1].ProteolysisProducts.Count());
            foreach (DisulfideBond bond in ok2[0].DisulfideBonds)
            {
                Assert.AreEqual(ok2[0].BaseSequence[bond.OneBasedBeginPosition - 1], 'C');
                Assert.AreEqual(ok2[0].BaseSequence[bond.OneBasedEndPosition - 1], 'C');
            }
            foreach (DisulfideBond bond in ok2[1].DisulfideBonds)
            {
                Assert.AreEqual(ok2[1].BaseSequence[bond.OneBasedBeginPosition - 1], 'C');
                Assert.AreEqual(ok2[1].BaseSequence[bond.OneBasedEndPosition - 1], 'C');
            }
        }

        [Test]
        public static void TestSlideDecoyXML()
        {
            var nice = new List<Modification>();
            var ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"disulfidetests.xml"), true, DecoyType.Slide, nice, false,
                new string[] { "exclude_me" }, out Dictionary<string, Modification> un);

            Assert.AreEqual("MALLVHFLPLLALLALWEPKPTQAFVKQHLCGPHLVEALYLVCGERGFFYTPKSRREVEDPQVEQLELGGSPGDLQTLALEVARQKRGIVDQCCTSICSLYQLENYCN", ok2[0].BaseSequence);
            Assert.AreEqual("MTKAEVLQLLAGLHLVHALYAVLGVRFFPYLPLSARWVPDPQQEFLKLHGCPPDLQELLLLVCREKGGFVTQKCRSECELPQVEQYENGCSNGLLYTSAIETACQDRI", ok2[1].BaseSequence);
            Assert.AreEqual(ok2[0].DisulfideBonds.Count(), ok2[1].DisulfideBonds.Count());
            Assert.AreEqual(ok2[0].ProteolysisProducts.Count(), ok2[1].ProteolysisProducts.Count());
            for (int i = 0; i < ok2[0].ProteolysisProducts.Count(); i++)
            {
                Assert.AreEqual(ok2[0].ProteolysisProducts.ToArray()[i].OneBasedBeginPosition, ok2[1].ProteolysisProducts.ToArray()[i].OneBasedBeginPosition);
                Assert.AreEqual(ok2[0].ProteolysisProducts.ToArray()[i].OneBasedEndPosition, ok2[1].ProteolysisProducts.ToArray()[i].OneBasedEndPosition);
            }
            foreach (DisulfideBond bond in ok2[0].DisulfideBonds)
            {
                Assert.AreEqual(ok2[0].BaseSequence[bond.OneBasedBeginPosition - 1], 'C');
                Assert.AreEqual(ok2[0].BaseSequence[bond.OneBasedEndPosition - 1], 'C');
            }
            foreach (DisulfideBond bond in ok2[1].DisulfideBonds)
            {
                Assert.AreEqual(ok2[1].BaseSequence[bond.OneBasedBeginPosition - 1], 'C');
                Assert.AreEqual(ok2[1].BaseSequence[bond.OneBasedEndPosition - 1], 'C');
            }
        }

        [Test]
        public static void TestReverseDecoyFasta()
        {
            List<Protein> prots = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"fasta.fasta"), true, DecoyType.Reverse, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, out var a);

            Assert.AreEqual("MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", prots[0].BaseSequence);
            Assert.AreEqual("MGGFGYLTRGQRKLAYVVDMATVTKRKAHETYTVADRIVNELFVKLVGRTEEYILGSIRKVGGRRALRRIAPKTIGQINDRLVKRHRKAGGKGLGKGGKGRGS", prots[1].BaseSequence);
        }

        [Test]
        public static void TestSlideDecoyFasta()
        {
            List<Protein> prots = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", @"fasta.fasta"), true, DecoyType.Slide, false,
                ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotGeneNameRegex,
                ProteinDbLoader.UniprotOrganismRegex, out var a);

            Assert.AreEqual("MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", prots[0].BaseSequence);
            Assert.AreEqual("MVRRRNAQGIGKGAGRKLRRSGGVGRGSKLLYKEGRKVHKKFLEDVIRGATTPTIHRKAKRVGAKDIVGAIKEQTRGLLGVGLGNFIYDTVGYRELAYRVTMT", prots[1].BaseSequence);
        }
    }
}