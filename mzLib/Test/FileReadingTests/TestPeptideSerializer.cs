using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics;
using Proteomics.ProteolyticDigestion;
using Proteomics;
using System.Diagnostics.CodeAnalysis;

namespace Test.FileReadingTests;

[ExcludeFromCodeCoverage]
public class TestPeptideSerializer
{
    // Peptide tests that were moved
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

        deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer>(), null);

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

        deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

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

        deserializedPeptide.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

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



    // Oligo Tests 

    //[Test]
    //public static void TestSerializationOligoFromString()
    //{
    //    // purpose of this test is to serialize/deserialize a OligoWithSetModifications and make sure the deserialized peptide
    //    // has the same properties as before it was serialized. This peptide is unmodified and generated from reading in a string
    //    string sequence = "PEPTIDE";
    //    OligoWithSetMods peptide = new OligoWithSetMods(sequence, new Dictionary<string, Modification>(), 0, null, null, 1, 7, 0);
    //    OligoWithSetMods deserializedPeptide = null;

    //    string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationOligoWithSetModsFromString");
    //    System.IO.Directory.CreateDirectory(dir);
    //    string path = System.IO.Path.Combine(dir, "myOligoIndex.ind");

    //    var messageTypes = typeof(OligoWithSetMods);
    //    var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

    //    using (var file = System.IO.File.Create(path))
    //    {
    //        ser.Serialize(file, peptide);
    //    }

    //    using (var file = System.IO.File.OpenRead(path))
    //    {
    //        deserializedPeptide = (OligoWithSetMods)ser.Deserialize(file);
    //    }

    //    deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer>(), null);

    //    // not asserting any protein properties - since the peptide was created from a sequence string it didn't have a protein to begin with

    //    Assert.That(peptide.Equals(deserializedPeptide));
    //    Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
    //    Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

    //    var products = new List<Product>();

    //    deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
    //    List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

    //    peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
    //    List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

    //    Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
    //}



    //[Test]
    //public static void TestSerializationOligoFromRna()
    //{
    //    // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
    //    // has the same properties as before it was serialized. This peptide is unmodified and generated from digesting a protein
    //    RNA protein = new("GUACUGAGUCUACUAGAUCA", "name", "Accession1", "", "");

    //    OligoWithSetMods peptide = protein.Digest(new RnaDigestionParams("RNase T1"), new List<Modification>(), new List<Modification>()).First();
    //    OligoWithSetMods deserializedPeptide = null;

    //    string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProtein");
    //    System.IO.Directory.CreateDirectory(dir);
    //    string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

    //    var messageTypes = typeof(OligoWithSetMods);
    //    var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

    //    using (var file = System.IO.File.Create(path))
    //    {
    //        ser.Serialize(file, peptide);
    //    }

    //    using (var file = System.IO.File.OpenRead(path))
    //    {
    //        deserializedPeptide = (OligoWithSetMods)ser.Deserialize(file);
    //    }
    //    deserializedPeptide.SetNonSerializedPeptideInfo(new Dictionary<string, Modification>(), new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

    //    Assert.That(peptide.DigestionParams.Equals(deserializedPeptide.DigestionParams));
    //    Assert.That(peptide.Equals(deserializedPeptide));
    //    Assert.That(deserializedPeptide.Parent.Name == peptide.Parent.Name);
    //    Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
    //    Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

    //    var products = new List<Product>();

    //    deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
    //    List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

    //    peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
    //    List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

    //    Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
    //}



    //[Test]
    //public static void TestSerializationOligoFromRnaWithMod()
    //{
    //    // purpose of this test is to serialize/deserialize a PeptideWithSetModifications and make sure the deserialized peptide
    //    // has the same properties as before it was serialized. This peptide is modified with a phosphorylation

    //    ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

    //    Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
    //    {
    //        { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
    //        { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
    //    };

    //    Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");

    //    Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 4, new List<Modification> { mod } } };

    //    RNA protein = new("GUACUGAGUCUACUAGAUCA", "name", "Accession1", "", "");

    //    OligoWithSetMods peptide = protein.Digest(new RnaDigestionParams("RNase T1"), new List<Modification>(), new List<Modification>()).First();
    //    OligoWithSetMods deserializedPeptide = null;

    //    string dir = System.IO.Path.Combine(TestContext.CurrentContext.TestDirectory, "TestSerializationPeptideFromProteinWithMod");
    //    System.IO.Directory.CreateDirectory(dir);
    //    string path = System.IO.Path.Combine(dir, "myPeptideIndex.ind");

    //    var messageTypes = typeof(OligoWithSetMods);
    //    var ser = new NetSerializer.Serializer(new List<Type> { messageTypes });

    //    using (var file = System.IO.File.Create(path))
    //    {
    //        ser.Serialize(file, peptide);
    //    }

    //    using (var file = System.IO.File.OpenRead(path))
    //    {
    //        deserializedPeptide = (OligoWithSetMods)ser.Deserialize(file);
    //    }

    //    Dictionary<string, Modification> stringToMod = new Dictionary<string, Modification> { { mods.Values.First().First().IdWithMotif, mods.Values.First().First() } };

    //    deserializedPeptide.SetNonSerializedPeptideInfo(stringToMod, new Dictionary<string, IBioPolymer> { { protein.Accession, protein } }, peptide.DigestionParams);

    //    Assert.That(peptide.Equals(deserializedPeptide));
    //    Assert.That(deserializedPeptide.Parent.Name == peptide.Parent.Name);
    //    Assert.That(deserializedPeptide.MonoisotopicMass == peptide.MonoisotopicMass);
    //    Assert.That(deserializedPeptide.SequenceWithChemicalFormulas == peptide.SequenceWithChemicalFormulas);

    //    var products = new List<Product>();

    //    deserializedPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
    //    List<double> deserializedPeptideFragments = products.Select(v => v.NeutralMass).ToList();

    //    peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
    //    List<double> peptideFragments = products.Select(v => v.NeutralMass).ToList();

    //    Assert.That(deserializedPeptideFragments.SequenceEqual(peptideFragments));
    //}
}