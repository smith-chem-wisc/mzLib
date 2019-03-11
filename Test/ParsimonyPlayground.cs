using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    public class ParsimonyPlayground
    {
        //[Test]
        //public static void MyStupidTest()
        //{
        //    List<Protein> prots = ProteinDbLoader.LoadProteinFasta(@"C:\Users\Michael Shortreed\Downloads\MUS_uniprot_canonical_170511.fasta", true, DecoyType.Reverse, false,
        //        ProteinDbLoader.UniprotAccessionRegex, ProteinDbLoader.UniprotFullNameRegex, ProteinDbLoader.UniprotNameRegex, ProteinDbLoader.UniprotGeneNameRegex,
        //        ProteinDbLoader.UniprotOrganismRegex, out var a);

        //    Dictionary<string, List<string>> peptideAccessionDictionary = new Dictionary<string, List<string>>();

        //    foreach (Protein protein in prots)
        //    {
        //        var digestedProtein = protein.Digest(new DigestionParams(maxMissedCleavages: 0), new List<Modification>(), new List<Modification>());

        //        foreach (PeptideWithSetModifications peptide in digestedProtein)
        //        {
        //            string peptideSequence = peptide.BaseSequence;

        //            if (peptideSequence.Length > 6 && peptideSequence.Length < 51)
        //            {
        //                if (peptideAccessionDictionary.ContainsKey(peptideSequence))
        //                {
        //                    if (!peptideAccessionDictionary[peptideSequence].Contains(peptide.Protein.Accession))
        //                    {
        //                        peptideAccessionDictionary[peptideSequence].Add(peptide.Protein.Accession);
        //                    }
        //                }
        //                else
        //                {
        //                    peptideAccessionDictionary.Add(peptideSequence, new List<string>() { peptide.Protein.Accession });
        //                }
        //            }
        //        }
        //    }

        //    List<string> proteinsWithAUniquePeptide = new List<string>();
        //    List<string> uniquePeptides = new List<string>();
        //    foreach (KeyValuePair<string, List<string>> item in peptideAccessionDictionary)
        //    {
        //        if (item.Value.Count == 1)
        //        {
        //            uniquePeptides.Add(item.Key);
        //            if (!proteinsWithAUniquePeptide.Contains(item.Value[0]))
        //            {
        //                proteinsWithAUniquePeptide.Add(item.Value[0]);
        //            }
        //        }
        //    }

        //    Dictionary<int, int> proteinsPerPeptide = new Dictionary<int, int>();

        //    foreach (KeyValuePair<string, List<string>> item in peptideAccessionDictionary)
        //    {
        //        int proteinCount = item.Value.Count;
        //        if (!proteinsPerPeptide.ContainsKey(proteinCount))
        //        {
        //            proteinsPerPeptide.Add(proteinCount, 1);
        //        }
        //        else
        //        {
        //            proteinsPerPeptide[proteinCount]++;
        //        }
        //    }

        //    Assert.IsTrue(false);
        //}

        //[Test]
        //public static void MyNewParsimonyClass()
        //{
        //    DigestionParams d = new DigestionParams();

        //    var proteinA = new Protein(sequence: "P", accession: "A");
        //    var proteinB = new Protein(sequence: "PE", accession: "B");

        //    ProteinList proteinList1 = new ProteinList(new List<Protein>() { proteinA });
        //    ProteinList proteinList2 = new ProteinList(new List<Protein>() { proteinB });
        //    ProteinList proteinList3 = new ProteinList(new List<Protein>() { proteinA });

        //    ParsimonySequence peptideOne = new ParsimonySequence(new PeptideWithSetModifications(sequence: "P", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
        //    ParsimonySequence peptideTwo = new ParsimonySequence(new PeptideWithSetModifications(sequence: "PE", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);

        //    ParsimonySequenceList peptideList1 = new ParsimonySequenceList(new List<ParsimonySequence>() { peptideOne });
        //    ParsimonySequenceList peptideList2 = new ParsimonySequenceList(new List<ParsimonySequence>() { peptideTwo });

        //    ProteinGroup pg1 = new ProteinGroup(proteinList1, peptideList1);
        //    ProteinGroup pg2 = new ProteinGroup(proteinList2, peptideList2);
        //    ProteinGroup pg3 = new ProteinGroup(proteinList3, peptideList1);

        //    ProteinGroup pg4 = pg1;

        //    Assert.IsTrue(pg1.Equals(pg3));
        //    Assert.IsFalse(pg1.Equals(pg2));

        //    pg3.AddToThisProteinGroup(pg1);
        //    Assert.AreEqual(1, pg3.ProteinList.Proteins.Count);
        //    Assert.AreEqual(1, pg3.ProteinList.AccessionNumbers.Count);
        //    Assert.AreEqual(1, pg3.PepWithSetModsList.ParsimonySequences.Count);
        //    Assert.AreEqual(1, pg3.PepWithSetModsList.FullSequenceList.Count);

        //    pg3.AddToThisProteinGroup(pg2);
        //    Assert.AreEqual(2, pg3.ProteinList.Proteins.Count);
        //    Assert.AreEqual(2, pg3.ProteinList.AccessionNumbers.Count);
        //    Assert.AreEqual(2, pg3.PepWithSetModsList.ParsimonySequences.Count);
        //    Assert.AreEqual(2, pg3.PepWithSetModsList.FullSequenceList.Count);

        //    Assert.IsFalse(pg3.Equals(pg1));

        //    Assert.IsTrue(pg1.ProteinList == pg4.ProteinList);
        //    Assert.IsTrue(pg1.ProteinList.Equals(pg4.ProteinList));
        //}

        //[Test]
        //public static void LetsTrySomeParsimony()
        //{
        //    var proteinK = new Protein(sequence: "PEPTKPEPKPEKP", accession: "Lysine");
        //    var proteinR = new Protein(sequence: "PEPTRPEPRPERP", accession: "Arginine");

        //    DigestionParams d = new DigestionParams();

        //    ParsimonySequence peptide_PEK = new ParsimonySequence(new PeptideWithSetModifications(sequence: "PEK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
        //    ParsimonySequence peptide_PER = new ParsimonySequence(new PeptideWithSetModifications(sequence: "PER", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
        //    ParsimonySequence peptide_PEPK = new ParsimonySequence(new PeptideWithSetModifications(sequence: "PEPK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
        //    ParsimonySequence peptide_PEPR = new ParsimonySequence(new PeptideWithSetModifications(sequence: "PEPR", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
        //    ParsimonySequence peptide_PEPTK = new ParsimonySequence(new PeptideWithSetModifications(sequence: "PEPTK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
        //    ParsimonySequence peptide_PEPTR = new ParsimonySequence(new PeptideWithSetModifications(sequence: "PEPTR", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);

        //    ParsimonySequence peptide_P = new ParsimonySequence(new PeptideWithSetModifications(sequence: "P", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);

        //    PeptideWithSetModifications test_peptide_PEK = new PeptideWithSetModifications(sequence: "PEK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_PER = new PeptideWithSetModifications(sequence: "PER", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_PEPK = new PeptideWithSetModifications(sequence: "PEPK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_PEPR = new PeptideWithSetModifications(sequence: "PEPR", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_PEPTK = new PeptideWithSetModifications(sequence: "PEPTK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_PEPTR = new PeptideWithSetModifications(sequence: "PEPTR", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);

        //    PeptideWithSetModifications test_peptide_P = new PeptideWithSetModifications(sequence: "P", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);

        //    List<Protein> theDatabase = new List<Protein>() { proteinK, proteinR };

        //    //Test 1
        //    List<PeptideWithSetModifications> observedPeptides = new List<PeptideWithSetModifications>() { test_peptide_PEPK, test_peptide_PEK };
        //    Parsimony P1 = new Parsimony();
        //    P1.AssignProteingroups(observedPeptides, theDatabase);
        //    Assert.AreEqual(1, P1.ProteinGroups.Count);
        //    Assert.AreEqual("Lysine", P1.ProteinGroups[0].ProteinList.Proteins[0].Accession);

        //    //Test 2
        //    observedPeptides = new List<PeptideWithSetModifications>() { test_peptide_PEK, test_peptide_PER };
        //    Parsimony P2 = new Parsimony();
        //    P2.AssignProteingroups(observedPeptides, theDatabase);
        //    Assert.AreEqual(2, P2.ProteinGroups.Count);

        //    List<string> testAccessions = new List<string>();
        //    foreach (ProteinGroup group in P2.ProteinGroups)
        //    {
        //        testAccessions.AddRange(group.ProteinList.AccessionNumbers);
        //    }
        //    CollectionAssert.AreEquivalent(new List<string>() { "Lysine", "Arginine" }, testAccessions);

        //    //Test 3
        //    observedPeptides = new List<PeptideWithSetModifications>() { test_peptide_P, test_peptide_PEK };
        //    Parsimony P3 = new Parsimony();
        //    P3.AssignProteingroups(observedPeptides, theDatabase);
        //    Assert.AreEqual(2, P3.ProteinGroups.Count);

        //    testAccessions = new List<string>();
        //    foreach (ProteinGroup group in P3.ProteinGroups)
        //    {
        //        testAccessions.AddRange(group.ProteinList.AccessionNumbers);
        //    }
        //    CollectionAssert.AreEquivalent(new List<string>() { "Lysine", "Arginine" }, testAccessions.Distinct());

        //    //Test 4
        //    observedPeptides = new List<PeptideWithSetModifications>() { test_peptide_PEPK, test_peptide_PEK, test_peptide_PEPTK };
        //    Parsimony P4 = new Parsimony();
        //    P4.AssignProteingroups(observedPeptides, theDatabase);
        //    Assert.AreEqual(1, P4.ProteinGroups.Count);
        //    Assert.AreEqual("Lysine", P4.ProteinGroups[0].ProteinList.Proteins[0].Accession);
        //}

        //[Test]
        //public static void BipartiteGraphExample()
        //{
        //    var proteinOne = new Protein(sequence: "DKEKHKIKLK", accession: "one");
        //    var proteinTwo = new Protein(sequence: "IK", accession: "two");
        //    var proteinThree = new Protein(sequence: "GK", accession: "three");
        //    var proteinFour = new Protein(sequence: "CKNK", accession: "four");
        //    var proteinFive = new Protein(sequence: "EKIK", accession: "five");
        //    var proteinSix = new Protein(sequence: "CKGK", accession: "six");
        //    var proteinSeven = new Protein(sequence: "AKFK", accession: "seven");
        //    var proteinEight = new Protein(sequence: "IK", accession: "eight");
        //    var proteinNine = new Protein(sequence: "CKNK", accession: "nine");

        //    List<Protein> theDatabase = new List<Protein>() { proteinOne, proteinTwo, proteinThree, proteinFour, proteinFive, proteinSix, proteinSeven, proteinEight, proteinNine };

        //    DigestionParams d = new DigestionParams();

        //    PeptideWithSetModifications test_peptide_One = new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Two = new PeptideWithSetModifications(sequence: "CK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Three = new PeptideWithSetModifications(sequence: "DK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Four = new PeptideWithSetModifications(sequence: "EK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Five = new PeptideWithSetModifications(sequence: "FK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Six = new PeptideWithSetModifications(sequence: "GK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Seven = new PeptideWithSetModifications(sequence: "HK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Eight = new PeptideWithSetModifications(sequence: "IK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Nine = new PeptideWithSetModifications(sequence: "LK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
        //    PeptideWithSetModifications test_peptide_Ten = new PeptideWithSetModifications(sequence: "NK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);

        //    List<PeptideWithSetModifications> observedPeptides = new List<PeptideWithSetModifications>() { test_peptide_One, test_peptide_Two, test_peptide_Three, test_peptide_Four, test_peptide_Five, test_peptide_Six, test_peptide_Seven, test_peptide_Eight, test_peptide_Nine, test_peptide_Ten };
        //    Parsimony P1 = new Parsimony();
        //    P1.AssignProteingroups(observedPeptides, theDatabase);
        //    Assert.AreEqual(1, P1.ProteinGroups.Count);
        //}

        [Test]
        public void TestParsimonySequenceEquality()
        {
            DigestionParams d = new DigestionParams();
            ParsimonySequence one = new ParsimonySequence(new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
            ParsimonySequence two = new ParsimonySequence(new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
            ParsimonySequence three = new ParsimonySequence(new PeptideWithSetModifications(sequence: "LK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);

            Assert.AreEqual(one, two);
            Assert.AreNotEqual(one, three);
        }

        [Test]
        public void TestParsimonySequenceAsDictionaryKey()
        {
            DigestionParams d = new DigestionParams();
            ParsimonySequence one = new ParsimonySequence(new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
            ParsimonySequence two = new ParsimonySequence(new PeptideWithSetModifications(sequence: "HK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
            ParsimonySequence three = new ParsimonySequence(new PeptideWithSetModifications(sequence: "LK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);

            Dictionary<ParsimonySequence, int> myDictionary = new Dictionary<ParsimonySequence, int>() { { one, 1 }, { two, 2 }, { three, 3 } };

            Assert.AreEqual(myDictionary[one], 1);
            Assert.AreEqual(myDictionary[two], 2);
            Assert.AreEqual(myDictionary[three], 3);

            Assert.AreNotEqual(myDictionary[one], myDictionary[two]);
            Assert.AreNotEqual(myDictionary[one], myDictionary[three]);
            Assert.AreNotEqual(myDictionary[two], myDictionary[three]);
            Assert.AreEqual(myDictionary[one], myDictionary[one]);

            ParsimonySequence four = new ParsimonySequence(new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);

            var ex = Assert.Throws<ArgumentException>(() => myDictionary.Add(four, 4));
            Assert.That(ex.Message, Is.EqualTo("An item with the same key has already been added."));

            Assert.AreEqual(1, myDictionary[four]);
        }

        [Test]
        public void TestProteinListEquality()
        {
            Protein proteinOne = new Protein(sequence: "DKEKHKIKLK", accession: "one");
            Protein proteinTwo = new Protein(sequence: "DKEKHKIKLK", accession: "one");
            Protein proteinThree = new Protein(sequence: "K", accession: "Two");
            ProteinList listOne = new ProteinList(new List<Protein>() { proteinOne });
            ProteinList listTwo = new ProteinList(new List<Protein>() { proteinTwo });
            ProteinList listThree = new ProteinList(new List<Protein>() { proteinThree });

            Assert.AreEqual(proteinOne, proteinTwo);
            Assert.AreNotEqual(proteinOne, proteinThree);

            Assert.AreEqual(listOne.Proteins, listTwo.Proteins);
            Assert.AreNotEqual(listOne.Proteins, listThree.Proteins);

            Assert.AreEqual(listOne.AccessionNumbers, listTwo.AccessionNumbers);
            Assert.AreNotEqual(listOne.AccessionNumbers, listThree.AccessionNumbers);

            Assert.AreEqual(listOne, listTwo);
            Assert.AreNotEqual(listOne, listThree);
        }

        [Test]
        public void TestPeptideListEquality()
        {
            DigestionParams d = new DigestionParams();
            ParsimonySequence psOne = new ParsimonySequence(new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
            ParsimonySequence psTwo = new ParsimonySequence(new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);
            ParsimonySequence psThree = new ParsimonySequence(new PeptideWithSetModifications(sequence: "LK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d), false);

            ParsimonySequenceList listOne = new ParsimonySequenceList(new List<ParsimonySequence>() { psOne });
            ParsimonySequenceList listTwo = new ParsimonySequenceList(new List<ParsimonySequence>() { psTwo });
            ParsimonySequenceList listThree = new ParsimonySequenceList(new List<ParsimonySequence>() { psThree });

            Assert.AreEqual(psOne, psTwo);
            Assert.AreNotEqual(psOne, psThree);

            Assert.AreEqual(listOne, listTwo);
            Assert.AreNotEqual(listOne, listThree);
        }

        [Test]
        public void TestProtease()
        {
            Protease pOne = new Protease("trypsin", CleavageSpecificity.Full, "1", "tOne", new List<DigestionMotif>());
            Protease pTwo = new Protease("trypsin", CleavageSpecificity.Full, "1", "tOne", new List<DigestionMotif>());
            Protease pThree = new Protease("argC", CleavageSpecificity.Full, "1", "aThree", new List<DigestionMotif>());

            Assert.AreEqual(pOne, pTwo);
            Assert.AreNotEqual(pOne, pThree);
        }

        [Test]
        public void TestPeptideWithSetMods()
        {
            DigestionParams d = new DigestionParams();
            PeptideWithSetModifications psOne = new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
            PeptideWithSetModifications psTwo = new PeptideWithSetModifications(sequence: "AK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);
            PeptideWithSetModifications psThree = new PeptideWithSetModifications(sequence: "NK", allKnownMods: new Dictionary<string, Modification>(), digestionParams: d);

            Assert.AreEqual(psOne, psTwo);
            Assert.AreNotEqual(psOne, psThree);
        }
    }
}