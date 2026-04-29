// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ModificationSerializationTests.cs) is part of Proteomics.
//
// Proteomics is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Proteomics is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Proteomics. If not, see <http://www.gnu.org/licenses/>.

using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.Omics.Modifications
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class ModificationSerializationTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestSerializationPeptideFromString()
        {
            // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
            // has the same properties as before it was serialized. This peptide is unmodified and generated from reading in a string
            string sequence = "PEPTIDE";
            PeptideWithSetModifications peptide = new PeptideWithSetModifications(sequence, new Dictionary<string, Modification>(), 0, null, null, 1, 7, 0);
            PeptideWithSetModifications deserializedPeptide = null;

            string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromString");
            System.IO.Directory.CreateDirectory(dir);
            string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

            var messageTypes = typeof(PeptideWithSetModifications);
            var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

            using (var file = System.IO.File.Create(path))
            {
                ser.Serialize(file, peptide);
            }

            using (var file = System.IO.File.OpenRead(path))
            {
                deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
            }

            deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, Protein>(), null);

            // not asserting any protein properties - since the peptide was created from a sequence string it didn't have a protein to begin with

            Assert.That(peptide.Equals(deserializedPeptide));
            Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
            Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

            var products = new List<Product>();

            deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

            Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        }

        [Test]
        public static void TestSerializationPeptideFromProtein()
        {
            // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
            // has the same properties as before it was serialized. This peptide is unmodified and generated from digesting a protein
            Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein");

            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).First();
            PeptideWithSetModifications deserializedPeptide = null;

            string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProtein");
            System.IO.Directory.CreateDirectory(dir);
            string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

            var messageTypes = typeof(PeptideWithSetModifications);
            var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

            using (var file = System.IO.File.Create(path))
            {
                ser.Serialize(file, peptide);
            }

            using (var file = System.IO.File.OpenRead(path))
            {
                deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
            }

            deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, Protein> { { protein.Accession, protein } }, peptide.DigestionParams);

            Assert.That(peptide.DigestionParams.Equals(deserializedPeptide.DigestionParams));
            Assert.That(peptide.Equals(deserializedPeptide));
            Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
            Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
            Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

            var products = new List<Product>();

            deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

            Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        }

        [Test]
        public static void TestSerializationPeptideFromProteinWithMod()
        {
            // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
            // has the same properties as before it was serialized. This peptide is modified with a phosphorylation

            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

            Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
            {
                { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
                { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
            };

            Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mod } } };

            Protein protein = new Protein("PEPTIDE", "Accession1", name: "MyProtein", oneBasedModifications: mods);

            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(v => v.AllModsOneIsNterminus.Count == 1).First();
            PeptideWithSetModifications deserializedPeptide = null;

            string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProteinWithMod");
            System.IO.Directory.CreateDirectory(dir);
            string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

            var messageTypes = typeof(PeptideWithSetModifications);
            var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

            using (var file = System.IO.File.Create(path))
            {
                ser.Serialize(file, peptide);
            }

            using (var file = System.IO.File.OpenRead(path))
            {
                deserializedPeptide = (PeptideWithSetModifications)ser.Deserialize(file);
            }

            Dictionary<string, Modification> stringToMod = new Dictionary<string, Modification> { { mods.Values.First().First().IdWithMotif, mods.Values.First().First() } };

            deserializedPeptide.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, Protein> { { protein.Accession, protein } }, peptide.DigestionParams);

            Assert.That(peptide.Equals(deserializedPeptide));
            Assert.That(deserializedPeptide.Protein.Name == peptide.Protein.Name);
            Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
            Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

            var products = new List<Product>();

            deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

            Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
        }
    }
}
